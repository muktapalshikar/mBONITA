
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
result <- get_reactome_analysis_result(analysis_id = "b74ae58e-1e3d-11ed-a0f6-5acc5b40ce9b", reactome_url = "https://gsa.reactome.org")

# save the result
saveRDS(result, file = "my_ReactomeGSA_result.rds")

# get the overview over all pathways
all_pathways <- pathways(result)
