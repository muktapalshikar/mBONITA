
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
result <- get_reactome_analysis_result(analysis_id = "85eb8188-1e3e-11ed-b0a0-8a3228bb25a5", reactome_url = "https://gsa.reactome.org")

# save the result
saveRDS(result, file = "my_ReactomeGSA_result.rds")

# get the overview over all pathways
all_pathways <- pathways(result)
