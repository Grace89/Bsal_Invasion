# This code was written by: G. V. DiRenzo
# If you have any questions, please send them to: grace.direnzo@gmail.com


# Code Objective:
  # To use a Latin hyper cube sampler to evenly sample distributions for parameter combinations used in simulations for REACTIVE management file


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
setwd("/Volumes/GVD/Yeti/Bsal/")



### Writing simulation parameters to tab delimited text file so that a BASH 
## script can read this file and set up jobs on cluster

library(lhs)




# 2. Set up information for sampler ---------------------------------------------------------------



# Number of samples in the parameter space
n_samples <- 10000


# Total number of parameters
n_params <- 27
  



# 3. Initialize sampler ---------------------------------------------------------------


lhs_raw <- randomLHS(n = n_samples, # Number of rows / samples
                     k = n_params)  # Number of columns or parameter variables

head(lhs_raw)

# Initialize data frame
param_combos <- setNames(data.frame(matrix(ncol = n_params, nrow = n_samples)),
                         c("n.sites",  
                           "n.pre.Bsal.pre.treat", "n.pre.Bsal.post.treat",
                           "n.post.Bsal.pre.treat", "n.post.Bsal.post.treat", 
                           "n.Bsal.sites",
                           "psi_bh", 
                           "c_H", "c_H_pro", "c_H_pro_react", 
                           "c_S_pro", "c_S_pro_react", 
                           "phi_Hb", "phi_Hb_pro",  "phi_Hb_pro_react", 
                           "delta_s_pro", "delta_s_pro_react",  
                           "delta_L_pro", "delta_L_pro_react",
                           "e_S_pro", "e_S_pro_react",
                           "e_L_pro", "e_L_pro_react",
                           "g_S_pro", "g_S_pro_react", 
                           "d_L_pro", "d_L_pro_react"))

head(param_combos)



# 4. Use LHS  ---------------------------------------------------------------




# Set variables that will change

param_combos$n.sites <- round(qunif(lhs_raw[, 1], min = 2, max = 1000))

# Make sure there aren't more Bsal infected sites than the total number of sites
for(i in 1:nrow(param_combos)){
  if(param_combos$n.sites[i] < 5){
    param_combos$n.Bsal.sites[i] <- round(qunif(lhs_raw[i, 2], min = 1, max = param_combos$n.sites[i]))
  } else {param_combos$n.Bsal.sites[i] <- round(qunif(lhs_raw[i, 2], min = 1, max = 5))}
}


param_combos$n.pre.Bsal.pre.treat <- round(qunif(lhs_raw[, 3], min =  2, max = 10))
param_combos$n.pre.Bsal.post.treat <- round(qunif(lhs_raw[, 4], min =  2, max = 10))

param_combos$n.post.Bsal.pre.treat <- round(qunif(lhs_raw[, 5], min =  2, max = 10))
param_combos$n.post.Bsal.post.treat <- round(qunif(lhs_raw[, 6], min =  2, max = 10))

param_combos$psi_bh <- qunif(lhs_raw[, 7], min =  0, max = 1)

param_combos$c_H <- qunif(lhs_raw[, 8], min =  0, max = 1)
param_combos$c_H_pro <- qunif(lhs_raw[, 9], min = param_combos$c_H, max = 1)# Has to be higher than c_H
param_combos$c_H_pro_react <- qunif(lhs_raw[, 10], min =  param_combos$c_H_pro, max = 1) # Has to be higher than c_H_pro

param_combos$c_S_pro <- qunif(lhs_raw[, 11], min =  0, max = 1)
param_combos$c_S_pro_react <- qunif(lhs_raw[, 12], min =  param_combos$c_S_pro, max = 1) # Has to be higher than c_S_pro

param_combos$phi_Hb <- qunif(lhs_raw[, 13], min =  0, max = 1)
param_combos$phi_Hb_pro <- qunif(lhs_raw[, 14], min =  param_combos$phi_Hb, max = 1) # Has to be higher than phi_Hb
param_combos$phi_Hb_pro_react <- qunif(lhs_raw[, 15], min =  param_combos$phi_Hb_pro, max = 1) # Has to be higher than phi_Hb_pro

param_combos$delta_s_pro       <- qunif(lhs_raw[, 16], min =  0, max = 1)
param_combos$delta_s_pro_react <- qunif(lhs_raw[, 17], min =  0, max = 1)

param_combos$delta_L_pro       <- qunif(lhs_raw[, 18], min =  0, max = param_combos$delta_s_pro)
param_combos$delta_L_pro_react <- qunif(lhs_raw[, 19], min =  0, max = param_combos$delta_s_pro_react) # Make sure that this value is less than delta_s_pro

param_combos$e_S_pro <- qunif(lhs_raw[, 20], min =  0, max = 1)
param_combos$e_S_pro_react <- qunif(lhs_raw[, 21], min =  param_combos$e_S_pro, max = 1) # Should be higher than without treatment

param_combos$e_L_pro <- qunif(lhs_raw[, 22], min =  0, max = 1)
param_combos$e_L_pro_react <- qunif(lhs_raw[, 23], min = param_combos$e_L_pro, max = 1) # Should be higher than without treatment

param_combos$g_S_pro <- qunif(lhs_raw[, 24], min =  0, max = 1)
param_combos$g_S_pro_react <- qunif(lhs_raw[, 25], min =  0, max = param_combos$g_S_pro) # Should be lower than without treatment

param_combos$d_L_pro <- qunif(lhs_raw[, 26], min =  0, max = 1)
param_combos$d_L_pro_react <- qunif(lhs_raw[, 27], min =  param_combos$d_L_pro, max = 1) # Should be higher than without treatment


head(param_combos)




# 5. Write file  ---------------------------------------------------------------



write.csv(param_combos, file = "./ParameterCombinations/parameter_combos_Bsal_pro_react.csv", row.names = FALSE)

# End Script
