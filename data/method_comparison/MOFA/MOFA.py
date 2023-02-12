
######################################################
## Template script to train a MOFA+ model in Python ##
######################################################

from mofapy2.run.entry_point import entry_point
import pandas as pd
import io
import requests # to download the online data

###############
## Load data ##
###############

data = pd.read_csv("concatMelt.csv", index_col=0)
data

###########################
## Initialise MOFA model ##
###########################


## (1) initialise the entry point ##
ent = entry_point()


## (2) Set data options ##
# - scale_groups: if groups have significantly different ranges, it is good practice to scale each group to unit variance
# - scale_views: if views have significantly different ranges, it is good practice to scale each view to unit variance
ent.set_data_options(
	scale_groups = True, 
	scale_views = True
)


# (3, option 2) Set data using a long data frame
ent.set_data_df(data)


## (4) Set model options ##
# - factors: number of factors. Default is K=10
# - likelihods: likelihoods per view (options are "gaussian","poisson","bernoulli"). 
# 		Default is None, and they are infered automatically
# - spikeslab_weights: use spike-slab sparsity prior in the weights? (recommended TRUE)
# - ard_factors: use ARD prior in the factors? (TRUE if using multiple groups)
# - ard_weights: use ARD prior in the weights? (TRUE if using multiple views)

# Simple (using default values)
ent.set_model_options()



## (5) Set training options ##
# - iter: number of iterations
# - convergence_mode: "fast", "medium", "slow". 
#		For exploration, the fast mode is good enough.
# - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# - freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
# - dropR2: minimum variance explained criteria to drop factors while training.
# 		Default is None, inactive factors are not dropped during training
# - gpu_mode: use GPU mode? this needs cupy installed and a functional GPU, see https://cupy.chainer.org/
# - verbose: verbose mode?
# - seed: random seed


# Advanced (using personalised values)
ent.set_train_options(
	convergence_mode = "fast",
	seed = 1,
	iter = 10
)


## (6, optional) Set stochastic inference options##
# Only recommended with very large sample size (>1e6) and when having access to GPUs
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.25 to 0.5)
# - start_stochastic: first iteration to apply stochastic inference (recommended > 5)

# Simple (using default values)
# ent.set_stochastic_options()

# Advanced (using personalised values)
# ent.set_stochastic_options(batch_size=0.5, learning_rate=0.75, forgetting_rate=0.5, start_stochastic=10)


####################################
## Build and train the MOFA model ##
####################################

# Build the model 
ent.build()

# Run the model
ent.run()

##################################################################
## (Optional) do dimensionality reduction from the MOFA factors ##
##################################################################

# ent.umap()
# ent.tsne()

####################
## Save the model ##
####################

outfile = "model.hdf5"

# - save_data: logical indicating whether to save the training data in the hdf5 file.
# this is useful for some downstream analysis in R, but it can take a lot of disk space.
ent.save(outfile, save_data=True)
