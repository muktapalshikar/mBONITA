#Query gene sets
c1a <- read.csv("GOID_to_Name_Mapping.csv",stringsAsFactors = FALSE)
c1g <- read.csv("Gene_Symbol_GOID_mapping.csv",stringsAsFactors = FALSE)
c1b <- read.csv("GOID.csv",stringsAsFactors = FALSE)





#Sorted gene list for CNV, MET and MRNA
fall <- read.csv("Top1000Genes.csv",stringsAsFactors = FALSE)
f1 <- fall$CNV
f2 <- fall$MET
f3 <- fall$MRNA

commongenes <- Reduce(intersect, list(f1, f2, f3)) #check f1,f2 and f3 should be the same, but differently ordered

#Steps of random walk
f1a <- f1
f2a <- f2
f3a <- f3
steps <- 1:length(f1a)




#Precomputed LogN factorial to speed up the computation
logfac1 <- read.csv("LogNFactorial.csv",stringsAsFactors = FALSE)
logfac1 <- logfac1$LogFactorial


#Precomputed list of log of binomial coefficients for each step 
N <- length(commongenes)
logchoose <- function(N= 0, n = 0){
  if(N == n | n == 0){
    return(0)
  }else{
    loga1 <- logfac1[N]-logfac1[n] -logfac1[N-n]
    return(loga1)
  }
  
}

pnext <- list()
for(i in steps){
  nmax <- min(i,N-i)
  tempest <- as.numeric()
  pchoose1 <- logchoose(N, i)
  pnv1 <- N-i
  tempest <- lapply(0:nmax, function(jj){
    pnextra <- (logchoose(pnv1,jj)+logchoose(i, i - jj)- pchoose1)
  })
  pnext[[i]] <- (tempest)
}



for(ii in 1:length(c1b$GOName)){
  
  c1gname <- c1b$GOName[ii]
  c1ga <- c1g$HGNCSymbol[c1g$GOID == c1b$GOID[ii]]
  g1 <- c1ga
  g1n <- c1gname
  target <- g1
  
  #Number of observed cancer genes when taking 3 features into consideration
  overlap <- numeric(length(steps))
  
  for(i in steps){
    combined <- unique(c(f1a[1:i],f2a[1:i],f3a[1:i]))
    overlap[i] <- sum(combined %in% target)
  }
  
  
  #Number of observed cancer genes when considering feautre 1 (F1 = CNV) only
  f1aonly <- numeric(length(steps))
  for(i in steps){
    f1aonly[i] <- sum(f1a[1:i] %in% target)
  }
  
  #Number of observed cancer genes when considering feature 2 (F2 = Methylation) only
  f2aonly <- numeric(length(steps))
  for(i in steps){
    f2aonly[i] <- sum(f2a[1:i] %in% target)
  }
  
  #Number of observed cancer genes when considering feature 3 (F3 = MRNA) only
  f3aonly <- numeric(length(steps))
  for(i in steps){
    f3aonly[i] <- sum(f3a[1:i] %in% target)
  }
  
  #Number of observed cancer genes when considering F1F2 only
  f1f2aonly <- numeric(length(steps))
  for(i in steps){
    f1f2aonly[i] <- sum(unique(c(f1a[1:i],f2a[1:i])) %in% target)
  }
  
  #Number of observed cancer genes when considering F1F3
  f1f3aonly <- numeric(length(steps))
  for(i in steps){
    f1f3aonly[i] <- sum(unique(c(f1a[1:i],f3a[1:i])) %in% target)
  }
  
  #Number of observed cancer genes when considering F2F3
  f2f3aonly <- numeric(length(steps))
  for(i in steps){
    f2f3aonly[i] <- sum(unique(c(f2a[1:i],f3a[1:i])) %in% target)
  }
  
  
  # Expected cancer gene contributed by F1 given F2 (F1|F2)
  N <- length(commongenes)
  K <- sum(unique(c(f2a)) %in% target)
  f1gf2 <- numeric(length(steps))
  for(i in steps){
    geneexp <- f2a[1:i]
    geneepi <- f1a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
        #print(estgene)
      })
      sum(unlist(smallerest1))
    })
    f1gf2[i] <- sum(unlist(tempest))+k1
  }
  
  # Expected cancer gene contributed by F1 given F3 (F1|F3)
  N <- length(commongenes)
  K <- sum(unique(c(f3a)) %in% target)
  f1gf3 <- numeric(length(steps))
  for(i in steps){
    geneexp <- f3a[1:i]
    geneepi <- f1a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
        })
      sum(unlist(smallerest1))
    })
    f1gf3[i] <- sum(unlist(tempest))+k1
  }
  
  # Expected cancer gene contributed by F2 given F1 (F2|F1)
  N <- length(commongenes)
  K <- sum(unique(c(f1a)) %in% target)
  f2gf1 <- numeric(length(steps))
  for(i in steps){
    geneexp <- f1a[1:i]
    geneepi <- f2a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
      })
      sum(unlist(smallerest1))
    })
    f2gf1[i] <- sum(unlist(tempest))+k1
  }
  
  
  
  # Expected cancer gene contributed by F2 given F3 (F2|F3)
  N <- length(commongenes)
  K <- sum(unique(c(f3a)) %in% target)
  f2gf3 <- numeric(length(steps))
  for(i in steps){
    geneexp <- f3a[1:i]
    geneepi <- f2a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
      })
      sum(unlist(smallerest1))
    })
    f2gf3[i] <- sum(unlist(tempest))+k1
  }
  
  # Expected cancer gene contributed by F3 given F1 (F3|F1)
  N <- length(commongenes)
  K <- sum(unique(c(f1a)) %in% target)
  f3gf1 <- numeric(length(steps))
  for(i in steps){
    geneexp <- f1a[1:i]
    geneepi <- f3a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
      })
      sum(unlist(smallerest1))
    })
    f3gf1[i] <- sum(unlist(tempest))+k1
  }
  
  # Expected cancer gene contributed by F3 given F2 (F3|F2)
  N <- length(commongenes)
  K <- sum(unique(c(f2a)) %in% target)
  f3gf2 <- numeric(length(steps))
  for(i in steps){
    geneexp <- f2a[1:i]
    geneepi <- f3a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
      })
      sum(unlist(smallerest1))
    })
    f3gf2[i] <- sum(unlist(tempest))+k1
  }
  
  # Expected cancer gene contributed by F1F2F3 given F1F2 (F1F2F3|F1F2)
  N <- length(commongenes)
  K <- sum(unique(c(f1a,f2a)) %in% target)
  f1f2 <- numeric(length(steps))
  for(i in steps){
    geneexp <- unique(c(f1a[1:i],f2a[1:i]))
    geneepi <- f3a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
      })
      sum(unlist(smallerest1))
    })
    f1f2[i] <- sum(unlist(tempest))+k1
  }
  
  # Expected cancer gene contributed by F1F2F3 given F1F3 (F1F2F3|F1F3)
  N <- length(commongenes)
  K <- sum(unique(c(f1a,f3a)) %in% target)
  f1f3 <- numeric(length(steps))
  for(i in steps){
    geneexp <- unique(c(f1a[1:i],f3a[1:i]))
    geneepi <- f2a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
      })
      sum(unlist(smallerest1))
    })
    f1f3[i] <- sum(unlist(tempest))+k1
  }
  
  # Expected cancer gene contributed by F1F2F3 given F2F3 (F1F2F3|F2F3)
  N <- length(commongenes)
  K <- sum(unique(c(f2a,f3a)) %in% target)
  f2f3 <- numeric(length(steps))
  for(i in steps){
    geneexp <- unique(c(f2a[1:i],f3a[1:i]))
    geneepi <- f1a[1:i]
    nextra <- sum(!(geneepi %in% geneexp)) 
    k1 <- sum(geneexp %in% target)
    nmax <- min(i,N-i)
    kmax <- min(nextra,K-k1)
    tempest <- lapply(0:nmax, function(jj){
      pnextra <- pnext[[i]][[jj+1]]
      pkext1 <- logchoose(N-i,jj)
      pkv1 <- K-k1
      pkv2 <- N-i-K+k1
      jj1 <- min(jj, kmax)
      smallerest1 <- lapply(0:jj1, function(kk){
        pkv3 <- jj-kk
        pkextra <- (logchoose(pkv1, kk)+logchoose(pkv2,pkv3)-pkext1)
        estgene <- kk*exp(pnextra+pkextra)
      })
      sum(unlist(smallerest1))
    })
    f2f3[i] <- sum(unlist(tempest))+k1
  }
  
  df1 <- data.frame(F1 <-f1aonly,F2 <- f2aonly,F3<-f3aonly,F1F2 <- f1f2aonly,F1F3 <- f1f3aonly,F2F3 <- f2f3aonly, F1gF2 <-f1gf2,F1gF3 <- f1gf3,F2gF1 <- f2gf1,F2gF3 <- f2gf3,F3gF1 <- f3gf1,F3gF2 <- f3gf2,F1gF2F3 <- f2f3, F2gF1F3 <- f1f3,F3gF1F2 <- f1f2, Overlap <- overlap,stringsAsFactors = FALSE)
  name1 <- gsub(" ", "", c1gname)
  dftitle <- paste(name1, "DF.csv",sep = "")
  write.csv(df1,dftitle,row.names = FALSE)
  
  #Vs control
  expest <- seq(0,max(overlap), by = (max(overlap))/(length(steps)-1))
  expstata1 <- wilcox.test(f1aonly, expest,alternative = "greater")$p.value
  expstata2 <- wilcox.test(f2aonly, expest,alternative = "greater")$p.value
  expstata3 <- wilcox.test(f3aonly, expest,alternative = "greater")$p.value
  expstata4 <- wilcox.test(f1f2aonly, expest,alternative = "greater")$p.value
  expstata5 <- wilcox.test(f2f3aonly, expest,alternative = "greater")$p.value
  expstata6 <- wilcox.test(f1f3aonly, expest,alternative = "greater")$p.value
  expstata7 <- wilcox.test(overlap, expest,alternative = "greater")$p.value
  
  #F1F2
  expstat1 <- wilcox.test(f1f2aonly, f2gf1,alternative = "greater")$p.value
  expstat2 <- wilcox.test(f1f2aonly, f1gf2,alternative = "greater")$p.value
  #F1F3
  expstat3 <- wilcox.test(f1f3aonly, f3gf1,alternative = "greater")$p.value
  expstat4 <- wilcox.test(f1f3aonly, f1gf3,alternative = "greater")$p.value
  #F2F3
  expstat5 <- wilcox.test(f2f3aonly, f3gf2,alternative = "greater")$p.value
  expstat6 <- wilcox.test(f2f3aonly, f2gf3,alternative = "greater")$p.value
  
  #F3|F1F2
  expstat7 <- wilcox.test(overlap, f1f2,alternative = "greater")$p.value
  #F2|F1F3
  expstat8 <- wilcox.test(overlap, f1f3,alternative = "greater")$p.value
  #F1|F2F3
  expstat9 <- wilcox.test(overlap, f2f3,alternative = "greater")$p.value
  
  
  title1 <- paste(name1, "PVal.txt",sep = "")
  vals <- as.character(c(expstata1,expstata2,expstata3,expstat1,expstat2,expstat3,expstat4,expstat5,expstat6,expstat7,expstat8,expstat9,expstata4,expstata5,expstata6,expstata7))
  writeLines(text = vals, title1)
  
  
  #Multivariate GSEA plot
  library(ggplot2)
  
  #CNV vs MET
  expest <- seq(0,max(overlap), by = (max(overlap))/(length(steps)-1))
  dat1 <- data.frame(
    Rank = (c(steps,steps,steps,steps,steps,steps)),
    CanGenes = (c( f1f2aonly,f1aonly,f2aonly,f2gf1,f1gf2,expest)),
    Feature = c(rep("Combined",length(steps)),rep("CNV", length(steps)),rep("MET", length(steps)),rep("MET|CNV", length(steps)),rep("CNV|MET", length(steps)),rep("Null", length(steps)))
  )
  dat1$Feature<- factor(dat1$Feature, levels = c("Combined","CNV", "MET","MET|CNV","CNV|MET","Null"))
  cnvmet1 <- ggplot(data=dat1, aes(x=Rank, y=CanGenes, group=Feature, color=Feature)) +
    geom_line() + theme(legend.text=element_text(size=8))+ylab("Feature Count") +
    ggtitle("CNV vs MET",subtitle = "A")+ scale_color_manual(values=c("red","yellow","magenta","green","blue","black"))+geom_line(size=1)+theme(plot.title = element_text(hjust = 0.5))
  
  # CNV&MET VS MRNA
  dat1 <- data.frame(
    Rank = (c(steps,steps,steps,steps)),
    CanGenes = (c( overlap,f1f2aonly,f3aonly,f1f2)),
    Feature = c(rep("Combined",length(steps)),rep("CNV&MET", length(steps)),rep("MRNA", length(steps)),rep("MRNA|CNV&MET", length(steps)))
  )
  dat1$Feature<- factor(dat1$Feature, levels = c("Combined","CNV&MET", "MRNA","MRNA|CNV&MET"))
  cnvmet2 <- ggplot(data=dat1, aes(x=Rank, y=CanGenes, group=Feature, color=Feature)) +
    geom_line() + theme(legend.text=element_text(size=8))+ylab("Feature Count") +
    ggtitle("CNV&MET vs MRNA",subtitle = "B")+ scale_color_manual(values=c("red","yellow","magenta","green"))+geom_line(size=1)+theme(plot.title = element_text(hjust = 0.5))
  
  
  #CNV VS MRNA
  dat1 <- data.frame(
    Rank = (c(steps,steps,steps,steps,steps,steps)),
    CanGenes = (c( f1f3aonly,f1aonly,f3aonly,f3gf1,f1gf3,expest)),
    Feature = c(rep("Combined",length(steps)),rep("CNV", length(steps)),rep("MRNA", length(steps)),rep("MRNA|CNV", length(steps)),rep("CNV|MRNA", length(steps)),rep("Null", length(steps)))
  )
  dat1$Feature<- factor(dat1$Feature, levels = c("Combined","CNV", "MRNA","MRNA|CNV","CNV|MRNA","Null"))
  cnvmrna1 <- ggplot(data=dat1, aes(x=Rank, y=CanGenes, group=Feature, color=Feature)) +
    geom_line() + theme(legend.text=element_text(size=8))+ylab("Feature Count") +
    ggtitle("CNV vs MRNA",subtitle = "E")+ scale_color_manual(values=c("red","yellow","magenta","green","blue","black"))+geom_line(size=1)+theme(plot.title = element_text(hjust = 0.5))
  
  #CNV&MRNA VS MET
  dat1 <- data.frame(
    Rank = (c(steps,steps,steps,steps)),
    CanGenes = (c( overlap,f1f3aonly,f2aonly,f1f3)),
    Feature = c(rep("Combined",length(steps)),rep("CNV&MRNA", length(steps)),rep("MET", length(steps)),rep("MET|CNV&MRNA", length(steps))))
  dat1$Feature<- factor(dat1$Feature, levels = c("Combined","CNV&MRNA", "MET","MET|CNV&MRNA"))
  cnvmrna2 <- ggplot(data=dat1, aes(x=Rank, y=CanGenes, group=Feature, color=Feature)) +
    geom_line() + theme(legend.text=element_text(size=8))+ylab("Feature Count") +
    ggtitle("CNV&MRNA vs MET",subtitle = "F")+ scale_color_manual(values=c("red","yellow","magenta","green"))+geom_line(size=1)+theme(plot.title = element_text(hjust = 0.5))
  
  #MET VS MRNA
  dat1 <- data.frame(
    Rank = (c(steps,steps,steps,steps,steps,steps)),
    CanGenes = (c( f2f3aonly,f2aonly,f3aonly,f3gf2,f2gf3,expest)),
    Feature = c(rep("Combined",length(steps)),rep("MET", length(steps)),rep("MRNA", length(steps)),rep("MRNA|MET", length(steps)),rep("MET|MRNA", length(steps)),rep("Null", length(steps)))
  )
  dat1$Feature<- factor(dat1$Feature, levels = c("Combined","MET", "MRNA","MRNA|MET","MET|MRNA","Null"))
  metmrna1 <- ggplot(data=dat1, aes(x=Rank, y=CanGenes, group=Feature, color=Feature)) +
    geom_line() + theme(legend.text=element_text(size=8))+ylab("Feature Count") +
    ggtitle("MET vs MRNA",subtitle = "C")+ scale_color_manual(values=c("red","yellow","magenta","green","blue","black"))+geom_line(size=1)+theme(plot.title = element_text(hjust = 0.5))
  
  #MET&MRNA VS CNV
  dat1 <- data.frame(
    Rank = (c(steps,steps,steps,steps)),
    CanGenes = (c( overlap,f2f3aonly,f1aonly,f2f3)),
    Feature = c(rep("Combined",length(steps)),rep("MET&MRNA", length(steps)),rep("CNV", length(steps)),rep("CNV|MET&MRNA", length(steps)))
  )
  dat1$Feature<- factor(dat1$Feature, levels = c("Combined","MET&MRNA", "CNV","CNV|MET&MRNA"))
  metmrna2 <- ggplot(data=dat1, aes(x=Rank, y=CanGenes, group=Feature, color=Feature)) +
    geom_line() + theme(legend.text=element_text(size=8))+ylab("Feature Count")+
    ggtitle("MET&MRNA vs CNV",subtitle = "D")+ scale_color_manual(values=c("red","yellow","magenta","green"))+geom_line(size=1)+theme(plot.title = element_text(hjust = 0.5))
  
  name1 <- gsub(" ", "", c1gname)
  
  #Corresponding plots
  title1 <- paste("CNVMET", name1, "1.pdf",sep = "")
  pdf(title1, width = 15,height = 15)
  print(cnvmet1)
  dev.off()
  
  title1 <- paste("CNVMET", name1, "2.pdf",sep = "")
  pdf(title1, width = 15,height = 15)
  print(cnvmet2)
  dev.off()
  
  
  title1 <- paste("CNVMRNA", name1, "1.pdf",sep = "")
  pdf(title1, width = 15,height = 15)
  print(cnvmrna1)
  dev.off()
  
  title1 <- paste("CNVMRNA", name1, "2.pdf",sep = "")
  pdf(title1, width = 15,height = 15)
  print(cnvmrna2)
  dev.off()
  
  title1 <- paste("METMRNA", name1, "1.pdf",sep = "")
  pdf(title1, width = 15,height = 15)
  print(metmrna1)
  dev.off()
  
  
  title1 <- paste("METMRNA", name1, "2.pdf",sep = "")
  pdf(title1, width = 15,height = 15)
  print(metmrna2)
  dev.off()
}

