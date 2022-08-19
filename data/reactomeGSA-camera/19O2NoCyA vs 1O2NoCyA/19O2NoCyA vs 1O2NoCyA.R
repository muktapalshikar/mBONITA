
# This script downloads your recent ReactomeGSA result
# into an R session
#
# Note: The result is only stored for a certain period of time
#       on the ReactomeGSA servers. Therefore, it is highly
#       recommended to store the result locally.

# install the ReactomeGSA package if not available
if (!require(ReactomeGSA)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("ReactomeGSA")
}

# load the package
library(ReactomeGSA)

# load the analysis result
result <- get_reactome_analysis_result(analysis_id = "a390f677-1e38-11ed-9c61-1ab703fb60dd", reactome_url = "https://gsa.reactome.org")

# save the result
saveRDS(result, file = "19O2NoCyA vs 1O2NoCyA.rds")

# get the overview over all pathways
all_pathways <- pathways(result)

require(limma)
proteomics <- read.csv("bonita_proteomics.txt")
