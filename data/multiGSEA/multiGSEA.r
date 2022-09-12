## ----bioconductor, eval=FALSE-------------------------------------------------
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# # The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
#
# BiocManager::install("multiGSEA")
#

## ----devtools, eval=FALSE-----------------------------------------------------
#
#  install.packages("devtools")
#  devtools::install_github("https://github.com/yigbt/multiGSEA")
#

setwd("C:\\Users\\mpalshikar\\Documents\\moBONITA\\data\\multiGSEA")

## ----load_mapping_library, results='hide', warning=FALSE, message=FALSE-------
require( "org.Hs.eg.db")

## ----load_multiGSEA_package, results='hide', message=FALSE, warning=FALSE-----
require( multiGSEA)
require( magrittr)

## ----load_omics_data----------------------------------------------------------

require(readxl)


runMultiGSEA <- function(de_mrna, de_prot, de_phospho){

    de_mrna <- de_mrna[,c("Identifier", "logFC", "P.Value", "adj.P.Val")]
    colnames(de_mrna) <- c('Symbol', 'logFC', 'pValue', 'adj.pValue')

    de_prot <- de_prot[,c("Identifier", "logFC", "P.Value", "adj.P.Val")]
    colnames(de_prot) <- c('Symbol', 'logFC', 'pValue', 'adj.pValue')

    de_phospho <- de_phospho[,c("Identifier", "logFC", "P.Value", "adj.P.Val")]
    colnames(de_phospho) <- c('Symbol', 'logFC', 'pValue', 'adj.pValue')

    ## ----rank_features, results='hide'--------------------------------------------

    # create data structure
    omics_data <- initOmicsDataStructure( layer = c("transcriptome",
                                                    "proteome",
                                                    "metabolome"))

    ## add transcriptome layer
    omics_data$transcriptome <- rankFeatures( de_mrna$logFC,
                                              de_mrna$pValue)
    names( omics_data$transcriptome) <- de_mrna$Symbol

    ## add proteome layer
    omics_data$proteome <- rankFeatures(de_prot$logFC, de_prot$pValue)
    names( omics_data$proteome) <- de_prot$Symbol

    ## add phosphoproteome layer
    omics_data$metabolome <- rankFeatures(de_phospho$logFC, de_phospho$pValue)
    names( omics_data$metabolome) <- de_phospho$Symbol

    ## ----omics_short--------------------------------------------------------------

    omics_short <- lapply( names( omics_data), function( name){
        head( omics_data[[name]])
    })
    names( omics_short) <- names( omics_data)
    omics_short

    ## ----calculate_enrichment, results='hide', message=FALSE, warning=FALSE-------

    databases <- c( "kegg")#, "reactome")
    layers <- names( omics_data)

    pathways <- getMultiOmicsFeatures( dbs = databases, layer = c('transcriptome', 'proteome'),
                                       returnTranscriptome = "SYMBOL",
                                       returnProteome = "SYMBOL",
                                       #returnMetabolome = "HMDB",
                                       useLocal = FALSE)
    pathways$metabolome <- pathways$proteome

    ## ----pathways_short-----------------------------------------------------------

    pathways_short <- lapply( names( pathways), function( name){
        head( pathways[[name]], 2)
    })
    names( pathways_short) <- names( pathways)
    pathways_short

    ## ----run_enrichment, results='hide', message=FALSE, warning=FALSE-------------

    # use the multiGSEA function to calculate the enrichment scores
    # for all omics layer at once.
    enrichment_scores <- multiGSEA( pathways, omics_data)


    ## ----combine_pvalues----------------------------------------------------------

    df <- extractPvalues( enrichmentScores = enrichment_scores,
                          pathwayNames = names( pathways[[1]]))

    df$combined_pval <- combinePvalues( df)
    df$combined_padj <- p.adjust( df$combined_pval, method = "BH")

    df <- cbind( data.frame( pathway = names( pathways[[1]])), df)


    ## view results
    #head(df)
    #na.omit(df[df$combined_padj < 0.05,])
    #hist(df$combined_padj)
    return(df)
}


## Contrast: 19O2NoCyA vs 1O2NoCyA
excel_sheets(paste(getwd(),"/19O2NoCyA vs 1O2NoCyA/19O2NoCyA vs 1O2NoCyA.xlsx", sep = ""))
de_mrna <- read_excel(paste(getwd(),"/19O2NoCyA vs 1O2NoCyA/19O2NoCyA vs 1O2NoCyA.xlsx", sep = ""), sheet = 2)
de_prot <- read_excel(paste(getwd(),"/19O2NoCyA vs 1O2NoCyA/19O2NoCyA vs 1O2NoCyA.xlsx", sep = ""), sheet = 3)
de_phospho <- read_excel(paste(getwd(),"/19O2NoCyA vs 1O2NoCyA/19O2NoCyA vs 1O2NoCyA.xlsx", sep = ""), sheet = 4)

contrast1 <- runMultiGSEA(de_mrna, de_prot, de_phospho)
contrast1$Contrast <- "X19O2NoCyA - X1O2NoCyA"

# Contrast: 1O2PlusCyA vs 19O2NoCyA
excel_sheets(paste(getwd(),"/1O2PlusCyA vs 19O2NoCyA/1O2PlusCyA vs 19O2NoCyA.xlsx", sep = ""))
de_mrna <- read_excel(paste(getwd(),"/1O2PlusCyA vs 19O2NoCyA/1O2PlusCyA vs 19O2NoCyA.xlsx", sep = ""), sheet = 2)
de_prot <- read_excel(paste(getwd(),"/1O2PlusCyA vs 19O2NoCyA/1O2PlusCyA vs 19O2NoCyA.xlsx", sep = ""), sheet = 3)
de_phospho <- read_excel(paste(getwd(),"/1O2PlusCyA vs 19O2NoCyA/1O2PlusCyA vs 19O2NoCyA.xlsx", sep = ""), sheet = 4)

contrast2 <- runMultiGSEA(de_mrna, de_prot, de_phospho)
contrast2$Contrast <- "X19O2NoCyA - X1O2PlusCyA"


# Contrast: 1O2NoCyA vs 1O2PlusCyA
excel_sheets(paste(getwd(),"/1O2NoCyA vs 1O2PlusCyA/1O2NoCyA vs 1O2PlusCyA.xlsx", sep = ""))
de_mrna <- read_excel(paste(getwd(),"/1O2NoCyA vs 1O2PlusCyA/1O2NoCyA vs 1O2PlusCyA.xlsx", sep = ""), sheet = 2)
de_prot <- read_excel(paste(getwd(),"/1O2NoCyA vs 1O2PlusCyA/1O2NoCyA vs 1O2PlusCyA.xlsx", sep = ""), sheet = 3)
de_phospho <- read_excel(paste(getwd(),"/1O2NoCyA vs 1O2PlusCyA/1O2NoCyA vs 1O2PlusCyA.xlsx", sep = ""), sheet = 4)

contrast3 <- runMultiGSEA(de_mrna, de_prot, de_phospho)
contrast3$Contrast <- "X1O2NoCyA - X1O2PlusCyA"

allResults <- rbind(contrast1, contrast2, contrast3)
allResults$negativelog10FDR <- (-1)*log10(allResults$combined_padj)
allResults$Pathway <- stringr::str_to_title(trimws(gsub("\\(KEGG\\)", "", allResults$pathway)))
allResults$Algorithm <- "multiGSEA"
colnames(allResults) <- gsub("metabolome", "phosphoproteome", colnames(allResults))

write.csv(allResults, 'multiGSEA_allResults.csv')
print(head(allResults))

require(ggplot2)
require(ggpubr)
allResults <- na.omit(allResults[allResults$combined_padj< 0.05,])

ggplot(allResults) + geom_point(aes(x=negativelog10FDR, y=Pathway, color=Contrast)) + facet_wrap(facets = 'Contrast') + theme_pubr() + xlab("-log10(BH-corrected p-value)") + xlim(0, max(allResults$negativelog10FDR)+0.5)
ggsave(filename = "multiGSEA_allResults.svg", width = 9, height = 3, units = 'in')
ggsave(filename = "multiGSEA_allResults.png", width = 9, height = 3, units = 'in', dpi=600)
