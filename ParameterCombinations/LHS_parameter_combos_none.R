# This code was written by: G. V. DiRenzo
# If you have any questions, please send them to: grace.direnzo@gmail.com


# Code Objective:
  # To use a Latin hyper cube sampler to evenly sample distributions for parameter combinations used in simulations


# These parameters are used for the Bsal invasion simulation


##################################
###### Table of Contents #########
##################################

# 1. Set working directory & Load libraries
# 2. Format the parameter estimates
# 3. Set up the information for the simulation 
# 4. Set up storage
# 5. Indicate MCMC settings 
# 6. Save the parameters=
# 7. Create empty arrays to save outputs 
# 8. Simulate the data
# 9. Bundle everything for the model run 
# 10. Run the models 
# 11. Save the model output
# 12. Process the output 
# 13. Save one file with all the relevant information


##################################
##################################
##################################



# 1. Set working directory & Load libraries ---------------------------------------------------------------


# Set working directory
setwd("/lustre/projects/ecosystems/pwrc/gdirenzo/Bsal/")
#setwd("/Volumes/GVD/Yeti/Bsal/")


libLocation <- c(.libPaths(),
                 "/home/gdirenzo/R/x86_64-redhat-linux-gnu-library/3.6",
                 "/opt/ohpc/pub/libs/gnu8/R/3.6.1/lib64/R/library"
)


### Writing simulation parameters to tab delimited text file so that a BASH 
## script can read this file and set up jobs on cluster

library(lhs, lib.loc = libLocation)




# 2. Set up information for sampler ---------------------------------------------------------------



# Number of samples in the parameter space
n_samples <- 10000


# Total number of parameters
n_params <- 13
  



# 3. Initialize sampler ---------------------------------------------------------------


lhs_raw <- randomLHS(n = n_samples, # Number of rows / samples
                     k = n_params)  # Number of columns or parameter variables

head(lhs_raw)

# Initialize data frame
param_combos <- setNames(data.frame(matrix(ncol = n_params, nrow = n_samples)),
                         c("n.sites", "n.post.yrs",  "n.pre.yrs",
                           "psi_bh", 
                           "c_H", "c_S", 
                           "phi_Hb", "delta_s", "delta_L",
                           "e_S", "e_L",
                           "g_S", "d_L"))

head(param_combos)



# 4. Use LHS  ---------------------------------------------------------------




# Set variables that will change

param_combos$n.sites <- round(qunif(lhs_raw[, 1], min = 2, max = 1000))

param_combos$n.post.yrs <- round(qunif(lhs_raw[, 2], min = 1, max = 10))

param_combos$n.pre.yrs <- round(qunif(lhs_raw[, 3], min =  1, max = 10))

param_combos$psi_bh <- qunif(lhs_raw[, 4], min =  0, max = 1)

param_combos$c_H <- qunif(lhs_raw[, 5], min =  0, max = 1)

param_combos$c_S <- qunif(lhs_raw[, 6], min =  0, max = 1)

param_combos$phi_Hb <- qunif(lhs_raw[, 7], min =  0, max = 1)

param_combos$delta_s <- qunif(lhs_raw[, 8], min =  0, max = 1)
param_combos$delta_L <- qunif(lhs_raw[, 9], min =  0, max = param_combos$delta_s)

param_combos$e_S <- qunif(lhs_raw[, 10], min =  0, max = 1)
param_combos$e_L <- qunif(lhs_raw[, 11], min =  0, max = 1)

param_combos$g_S <- qunif(lhs_raw[, 12], min =  0, max = 1)
param_combos$d_L <- qunif(lhs_raw[, 13], min =  0, max = 1)


head(param_combos)




# 5. Write file  ---------------------------------------------------------------



write.csv(param_combos, file = "./ParameterCombinations/parameter_combos_Bsal_none.csv", row.names = FALSE)

# End Script