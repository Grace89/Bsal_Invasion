# This code was written by: G. V. DiRenzo
# If you have any questions, please send them to: grace.direnzo@gmail.com


# Code Objective:
  # To simulate Bsal invasion dynamics when there is only PROACTIVE management

# Output
  # Proportion of population in each of 6 states over time


# There are 6 states:
# 1 = bh = No Bsal, no host
# 2 = bH = No Bsal, host is present
# 3 = sh = Small Bsal amount, no host
# 4 = sH = Small Bsal amount, host is present
# 5 = Lh = Large Bsal amount, no host
# 6 = LH = Large Bsal amount, host is present


##################################
###### Table of Contents #########
##################################

# 1. Set working directory & Load libraries
# 2. HPC index #?# What is this?
# 3. Format the parameter estimates
# 4. Define the site/survey conditions
# 5. Define the initial and transition matrix 
# 6. Simulate data 
# 7. Summarize output
# 8. Save output


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


library("reshape2", lib.loc = libLocation)





# 2. HPC index ---------------------------------------------------------------


# temp
# niter <- 1

# Read in command line arguments
args<-commandArgs(trailingOnly=TRUE)
cat("the arguments are:")
print(args)

# Iteration number
niter <- as.numeric(args[1])




# 3. Format the parameter estimates ---------------------------------------------------------------



# Read in parameter combos
param_combos <- read.csv("./ParameterCombinations/parameter_combos_Bsal_pro.csv")


# Before Bsal invasion, the sites can only be in 1 of 2 states
  # Either bh or bH
psi_bh <- param_combos$psi_bh[niter]
psi_bH <- 1 - psi_bh

# Before Bsal invasion & before treatment - Probability of host colonizing a site
c_H    <- param_combos$c_H[niter]

# After proactive managemnet & before Bsal invasion - Probability of host colonizing a site
c_H_pro    <- param_combos$c_H_pro[niter]

# After proactive managemnet & Bsal invasion -  Probability of Bsal colonizing a site
c_S_pro    <- param_combos$c_S_pro[niter]

# Before Bsal invasion & before treatment - Host survival probability with NO Bsal & no proactive management
phi_Hb <- param_combos$phi_Hb[niter]

# After proactive managemnet - Host survival probability with NO Bsal
phi_Hb_pro <- param_combos$phi_Hb_pro[niter]

# After proactive managemnet & Bsal arrival - Multiplier for sites with Bsal present (small and large)
delta_s_pro <- param_combos$delta_s_pro[niter]
delta_L_pro <- param_combos$delta_L_pro[niter]

phi_Hs_pro <- phi_Hb_pro * delta_s_pro
phi_HL_pro <- phi_Hb_pro * delta_L_pro

# After proactive managemnet - Bsal extinction probability at a site where it is small
e_S_pro <- param_combos$e_S_pro[niter]

# After proactive managemnet - Bsal extinction probability at a site where it is large
e_L_pro <- param_combos$e_L_pro[niter]

# After proactive managemnet - Probability Bsal grows from small to large
g_S_pro <- param_combos$g_S_pro[niter]

# After proactive managemnet - Probability Bsal declines from large to small
d_L_pro <- param_combos$d_L_pro[niter]







# 4. Define the site/survey conditions ---------------------------------------------------------------



# Number of sites 
n.sites <- param_combos$n.sites[niter]

# Number of sites infected with Bsal when it arrives
n.Bsal.sites <- param_combos$n.Bsal.sites[niter]

# Number of states
n.states <- 6

# Time after Bsal arrivel but pre-treatment
n.pre.Bsal.pre.treat <- param_combos$n.pre.Bsal.pre.treat[niter]

# Time after Bsal arrivel and after pre-treatment
n.pre.Bsal.post.treat <- param_combos$n.pre.Bsal.post.treat[niter]

# Number of timesteps AFTER Bsal arrival
n.post.yrs <- param_combos$n.post.yrs[niter]

# Number of TOTAL timesteps
n.occasions <- n.pre.Bsal.pre.treat + n.pre.Bsal.post.treat + n.post.yrs




# 5. Define the initial and transition matrix ---------------------------------------------------------------



# Define the population
  # Initial conditions at timestep = 1
init_conditions <- matrix(NA, nrow = n.states, ncol = 1)

init_conditions[1] <- psi_bh
init_conditions[2] <- psi_bH

init_conditions[3] <- 0
init_conditions[4] <- 0

init_conditions[5] <- 0
init_conditions[6] <- 0


# Transition matrix after time step 1
  # Sites can be either bh or bH
  # Pre-Bsal arrival
trans_mat <- matrix(c(
  (1-c_H), c_H,
  (1-phi_Hb), phi_Hb
), nrow = 2, ncol = 2, byrow = TRUE)

# Add column & row names
colnames(trans_mat) <- rownames(trans_mat) <- c("bh", "bH")

# Make sure that all rows sum to 1
# rowSums(trans_mat)


# New matrix pre-Bsal invasion & after proactive treatment application
trans_mat_pro <- matrix(c((1-c_H_pro), c_H_pro,
                         (1-phi_Hb_pro), phi_Hb_pro), nrow = 2, ncol = 2, byrow = TRUE)

# Add column & row names
colnames(trans_mat_pro) <- rownames(trans_mat_pro) <- c("bh", "bH")



# After Bsal invasion
  # A site moves between these states with probabilities defined in the probability matrix
  # Defining the transition matrix
Bsal_trans_mat_pro <- matrix(c(
  (1-c_S_pro)*(1-c_H_pro),  (1-c_S_pro)*c_H_pro, c_S_pro*(1-c_H_pro), c_S_pro*c_H_pro, 0, 0,
  (1-c_S_pro)*(1-phi_Hb_pro), (1-c_S_pro)*phi_Hb_pro, c_S_pro*(1-phi_Hb_pro), c_S_pro*phi_Hb_pro, 0, 0,
  e_S_pro*(1-c_H_pro), e_S_pro*c_H_pro, (1-e_S_pro)*(1-g_S_pro)*(1-c_H_pro), (1-e_S_pro)*(1-g_S_pro)*c_H_pro, (1-e_S_pro)*g_S_pro*(1-c_H_pro), (1-e_S_pro)*g_S_pro*c_H_pro,
  e_S_pro*(1-phi_Hs_pro), e_S_pro*phi_Hs_pro , (1-e_S_pro)*(1-g_S_pro)*(1-phi_Hs_pro), (1-e_S_pro)*(1-g_S_pro)*phi_Hs_pro, (1-e_S_pro)*g_S_pro*(1-phi_Hs_pro), (1-e_S_pro)*g_S_pro*phi_Hs_pro,
  e_L_pro*(1-c_H_pro),    e_L_pro*c_H_pro, (1-e_L_pro)*d_L_pro*c_H_pro, (1-e_L_pro)*d_L_pro*(1-c_H_pro), (1-e_L_pro)*(1-d_L_pro)*c_H_pro, (1-e_L_pro)*(1-d_L_pro)*(1-c_H_pro),
  e_L_pro*(1-phi_HL_pro), e_L_pro*phi_HL_pro, (1-e_L_pro)*d_L_pro*(1-phi_HL_pro), (1-e_L_pro)*d_L_pro*phi_HL_pro, (1-e_L_pro)*(1-d_L_pro)*(1-phi_HL_pro), (1-e_L_pro)*(1-d_L_pro)*phi_HL_pro
), 
nrow = n.states, ncol = n.states, byrow = TRUE)


# Add column & row names
colnames(Bsal_trans_mat_pro) <- rownames(Bsal_trans_mat_pro) <- c("bh", "bH", "sh", "sH", "Lh", "LH")

# Make sure that all the rows sum to 1
rowSums(Bsal_trans_mat_pro)



# 6. Simulate data ---------------------------------------------------------------



# PROACTIVE management ---------------------------------------------------------------





# Create an empty vector to store probabilities
z <- array(NA, dim = c(n.sites, n.occasions))

# Randomly pick the sites where Bsal invades
infected.sites <- sample(1:n.sites, n.Bsal.sites)

# Use a for loop to determine what state each site is in
for(i in 1:n.sites){
  
  # First season
  z[i, 1] <- which(rmultinom(1, 1, init_conditions) == 1)
  
  # Before Bsal invasion & before proactive treatment
  for(k in 2:n.pre.Bsal.pre.treat){
    z[i,k] <- which(rmultinom(1, 1, trans_mat[z[i,k-1], ]) == 1)
  }
  
  ######### Treatment applied
  
  # Before Bsal invasion & after proactive treatment
  for(k in (n.pre.Bsal.pre.treat+1):n.pre.Bsal.post.treat){
    z[i,k] <- which(rmultinom(1, 1, trans_mat_pro[z[i,k-1], ]) == 1)
  }
  
}

######## Bsal introduction

# Determine if the sites that will get Bsal introduced have a host present or absent
for(i in 1:length(infected.sites)){
  
  # If a host is absent (bh = state 1) - the site will convert to Lh (state 5)
  if(z[infected.sites[i], n.pre.Bsal.post.treat] == 1){z[infected.sites[i], n.pre.Bsal.post.treat] <- 5}
  
  # If a host is present (bH = state 2) - the site will convert to LH (state 6)
  if(z[infected.sites[i], n.pre.Bsal.post.treat] == 2){z[infected.sites[i], n.pre.Bsal.post.treat] <- 6}
}

# Then - after Bsal arrival

for(i in 1:n.sites){  

  for(k in (n.pre.Bsal.post.treat+1):n.occasions){
    
    z[i,k] <- which(rmultinom(1, 1, Bsal_trans_mat_pro[z[i,k-1],]) ==1)
  
  }

}





# 7. Summarize output ---------------------------------------------------------------




# Create an empty matrix to summarize the proportion of sites in each state at each time step
Bsal_sum <- array(NA, dim = c(n.states, ncol(z)))


# Loop through each season and state to determine the proportion of sites
for(i in 1:n.occasions){
  
  for(k in 1:n.states){
    
    Bsal_sum[k,i] <- length(which(z[,i] == k)) / n.sites
    
  }
  
}


# Add rownames
rownames(Bsal_sum) <- c("bh", "bH", "sh", "sH", "Lh", "LH")

# Turn into long format
Bsal_sum.long <- melt(Bsal_sum)

# Add column names
colnames(Bsal_sum.long) <- c("state", "n.occasion", "proportion")





# 8. Create a summary of initial conditions, parameters, & outcome ---------------------------------------------------------------



# Subset the data from the last occasion
last.occ <- matrix(Bsal_sum.long$proportion[which(Bsal_sum.long$n.occasion == n.occasions)],
                   nrow = 1, byrow = TRUE)


# Add column names
colnames(last.occ) <-  c("bh", "bH", "sh", "sH", "Lh", "LH")


# Column bind the parameter combinations and last occ info together
params_dat <- cbind(param_combos[niter,], last.occ)





# 9. Save output ---------------------------------------------------------------



# Save the entire trajectory
save(Bsal_sum.long, file = paste0("./SimOutput/Trajectory/Bsal_sim_output_pro_", niter, ".rda"))



# Save the entire trajectory
save(params_dat, file = paste0("./SimOutput/Summary/Bsal_params_lastOcc_pro_", niter, ".rda"))


# End script
