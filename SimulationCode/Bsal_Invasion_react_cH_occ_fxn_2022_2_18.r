########################################
########################################
# This code was written by: G. V. DiRenzo & E. H. Grant
# If you have any questions, please email: gdirenzo@umass.edu, ehgrant@usgs.gov
########################################
########################################




########################################
####### Code Objective #################
########################################


# To simulate Bsal invasion host-pathogen dynamics

# In this scenario, we are evaluating the effectiveness of reactive actions to Bsal invasion


########################################
####### Code Output ####################
########################################



# This code generates the following files:
  # a summary figure of % of hosts remaining (total) over time
  # a summary figure of % of hosts in each group over time by scenario



########################################
########## Notes #######################
########################################




########################################
######## Table of Contents #############
########################################


# 1. Set working directory & Load libraries
# 2. Define the site/survey conditions
# 3. Format the parameter estimates
# 4. Define the function to run the simulation
# 5. Run the different scenarios
# 6. Create a figure with % of hosts in each group over time by scenario
# 7. Create a summary of % of hosts remaining in each group


########################################
########################################
########################################




# 1. Set working directory & Load libraries ---------------------------------------------------------------



# Remove objects in your working space
rm(list = ls(all.names = TRUE))
options(digits=3)

# Set working directory
# Evan
#setwd("C:\\Users\\ehgrant\\Documents\\Diseases\\WNS_Bsal Postdoc\\BSAL\\Transition model")

# Grace
#setwd("~/Dropbox/USGS/Bsal elicitation/")

# Molly
setwd("~/My_FILES/1_USGS_BsalDS_Postdoc/Patux_Transition_Mod")

#load libraries
library("reshape2")
library("tidyverse")
library("mc2d")


# 2. Define the site/survey conditions ---------------------------------------------------------------

# Number of simulations
n.sim <- 1

# Number of sites 
n.sites <- 10 #param_combos$n.sites[niter]

# Number of sites where Bsal arrives
n.Bsal.sites <- 1 #param_combos$n.Bsal.sites[niter]

# Number of states
n.states <- 6

# Timesteps before Bsal arrival
n.pre.yrs <- 5 #param_combos$n.pre.yrs[niter]




# 3. Format the parameter estimates ---------------------------------------------------------------


# There are 6 states each patch can be in: (#?# patch = site?-- Yes)
  # 1 = bh = No Bsal, no host
  # 2 = bH = No Bsal, host is present
  # 3 = sh = Small Bsal prevalence, no host
  # 4 = sH = Small Bsal prevalence, host is present
  # 5 = Lh = Large Bsal prevalence, no host
  # 6 = LH = Large Bsal prevalence, host is present


# Read in parameter combos
# param_combos <- read.csv("./parameter_combos_Bsal_react.csv") 

# Before Bsal invasion, the sites can only be in 1 of 2 states
# Either bh or bH
# psi_bh = Probability a site is not occupied by a host (state = bh) before Bsal arrival
psi_bh <- rpert(n.sim, min = 0, max = 0.40, mode = 0.2) 
    ##rpert random number from the PERT distribution - frequently used (with the triangular distribution) to translate expert estimates of the min, max and mode of a random variable in a smooth parametric distribution. #?# Why mode and not mean -- look into?

# psi_bH = Probability a site is occupied by a host (state = bH) before Bsal arrival
psi_bH <- 1 - psi_bh  #?# is there a difference between defining "not occ" with rpert vs occ  Can't do rpert again because has to equal 1; doing bh first is based on questions?

# Before any management & Before Bsal invasion- Probability of host colonizing a site
c_H    <- rpert(n.sim, min = 0.10, max = 0.50, mode = 0.2) 

# Before any management & after Bsal invasion - Probability of Bsal colonizing a site
c_S    <- rpert(n.sim, min = 0.60, max = 0.95, mode = 0.8) 

# After reactive management & after Bsal invasion - Probability of Bsal colonizing a site
c_S_react  <- rpert(n.sim, min = 0.01, max = 0.3, mode = 0.1) 

# After reactive management & after Bsal invasion - Probability of host colonizing a site
#c_H_react <- c_H

c_H_alpha <- 0 # intercept

c_H_beta <- 0.25 # slope

# All time before any management (includes time before and following Bsal arrival) - Host survival probability with NO Bsal
phi_Hb <- rpert(n.sim, min = 0.80, max = 0.98, mode = 0.95) #?# is this survival probability of an individual or the site/population?  the probability that the site remains occupied

# After reactive management & after Bsal invasion - Host survival probability with NO Bsal
phi_Hb_react <- rpert(n.sim, min = 0.80, max = 0.98, mode = 0.95) 

## Before any management & after Bsal invasion - Multiplier for sites with Bsal present (small and large)
phi_Hs <- rpert(n.sim, min = 0.01, max = 0.20, mode = 0.08) 
phi_HL <- rpert(n.sim, min = 0.001, max = 0.15, mode = 0.05) 

## After reactive management & after Bsal invasion - Multiplier for sites with Bsal present (small and large)
  # Management does not affect the baseline susceptibility  #?# reactive management doesn't but proactive could right? -- Yes
phi_Hs_react <- phi_Hs
phi_HL_react <- phi_HL

# Before any management & after Bsal invasion - Bsal extinction probability at a site where it is small
e_S <- rpert(n.sim, min = 0.01, max = 0.1, mode = 0.05) 

# After reactive management & after Bsal invasion - Bsal extinction probability at a site where it is small
e_S_react <- rpert(n.sim, min = 0.60, max = 0.80, mode = 0.7) 

# Before any management & after Bsal invasion - Bsal extinction probability at a site where it is large
e_L <- rpert(n.sim, min = 0.001, max = 0.1, mode = 0.01) 

# After reactive management & after Bsal invasion - Bsal extinction probability at a site where it is large
e_L_react <- rpert(n.sim, min = 0.05, max = 0.35, mode = 0.15) 

# Before any management & after Bsal invasion - Probability Bsal grows from small to large
g_S <- rpert(n.sim, min = 0.60, max = 0.95, mode = 0.8) 

# After reactive management & after Bsal invasion - Probability Bsal grows from small to large
g_S_react <- rpert(n.sim, min = 0.01, max = 0.3, mode = 0.1) 

# Before any management & after Bsal invasion - Probability Bsal declines from large to small
d_L <- rpert(n.sim, min = 0.01, max = 0.1, mode = 0.06) 

# After reactive management & after Bsal invasion - Probability Bsal declines from large to small
d_L_react <- rpert(n.sim, min = 0.1, max = 0.35, mode = 0.2) 




# 4. Define the function to run the simulation ---------------------------------------------------------------


Bsal.sim.react <- function(
  
  # Number of timesteps post Bsal arrival & before treatment
  n.yrs.post.Bsal.pre.treat, 
  
  # Number of timesteps post Bsal arrival & after treatment     
  n.post.Bsal.post.treat, 
  
  # total number of timesteps
  n.occasions){ 
  
  # Define the population
  # Initial conditions at timestep = 1
    # Here- we use information about host occupancy probability (psi_bh and psi_bH) - no hosts are in states 3 - 6 (which correspond to states after Bsal arrival)
  init_conditions <- matrix(NA, nrow = n.states, ncol = n.sim)
  
  # All columns
  init_conditions[1, ] <- psi_bh  # row 1  # prob a site is unoccupied
  init_conditions[2, ] <- psi_bH  # row 2  # prob a site is occupied
  
  init_conditions[3, ] <- 0       # row 3
  init_conditions[4, ] <- 0       # row 4
  
  init_conditions[5, ] <- 0       # row 5
  init_conditions[6, ] <- 0       # row 6
  
  
  # Transition matrix after time step 1
  # Sites can be either bh or bH
    # No Bsal
    # Site is occupied by the host or not
  # Sites can be colonized by the host (c_H)
  # A host can leave a site - making it go "extinct" - 1 - phi_Hb
  # Or a site can remain occupied phi_Hb 
  
  trans_mat <- array(NA, dim = c(2, 2,n.occasions, n.sim)) # now has 4 dimensions because cH is a function of occ at a given time
  
  
  
  # Before any treatment & After Bsal invasion
  # A site moves between these states with probabilities defined in the probability matrix
  # Defining the transition matrix
  Bsal_trans_mat <- array(NA, dim = c(n.states, n.states, n.occasions, n.sim))

  
  # Check that all the rows sum to 1
  # rowSums(Bsal_trans_mat[, , 1])
  
  
  # After reactive treatment & Bsal arrival
  # A site moves between these states with probabilities defined in a new probability matrix
  # Defining the transition matrix
  Bsal_trans_mat_react <- array(NA, dim = c(n.states, n.states,n.occasions, n.sim))
  
  for(i in 1:n.sim){
    Bsal_trans_mat_react[,,i] <- matrix(c(
      
      # Row 1
      (1-c_S_react[i])*(1-c_H_react[i]),  # Col 1  bh -> bh
      (1-c_S_react[i])*c_H_react[i],      # Col 2  bh -> bH
      c_S_react[i]*(1-c_H_react[i]),      # Col 3  bh -> sh
      c_S_react[i]*c_H_react[i],          # Col 4  bh -> sH
      0,                                  # Col 5  bh -> Lh
      0,                                  # Col 6  bh -> LH
      
      # Row 2
      (1-c_S_react[i])*(1-phi_Hb_react[i]), # Col 1  bH -> bh #?# are we asking experts for both survival prob if Bsal is there and its not?? Yes based on l151-157  No but with get it from the 1 - its there.
      (1-c_S_react[i])*phi_Hb_react[i],     # Col 2  bH -> bH
      c_S_react[i]*(1-phi_Hb_react[i]),     # Col 3  bH -> sh
      c_S_react[i]*phi_Hb_react[i],         # Col 4  bH -> sH
      0,                                    # Col 5  bH -> Lh
      0,                                    # Col 6  bH -> LH
      
      # Row 3
      e_S_react[i]*(1-c_H_react[i]),                      # Col 1 sh -> bh
      e_S_react[i]*c_H_react[i],                          # Col 2 sh -> bH
      (1-e_S_react[i])*(1-g_S_react[i])*(1-c_H_react[i]), # Col 3 sh -> sh
      (1-e_S_react[i])*(1-g_S_react[i])*c_H_react[i],     # Col 4 sh -> sH
      (1-e_S_react[i])*g_S_react[i]*(1-c_H_react[i]),     # Col 5 sh -> Lh **** Does not match trans matrix written--> MCB updated PP 1/18/22
      (1-e_S_react[i])*g_S_react[i]*c_H_react[i],         # Col 6 sh -> LH **** Does not match trans matrix written--> MCB updated PP 1/18/22
      
      # Row 4
      e_S_react[i]*(1-phi_Hs_react[i]),                     # Col 1  sH -> bh
      e_S_react[i]*phi_Hs_react[i],                         # Col 2  sH -> bH
      (1-e_S_react[i])*(1-g_S_react[i])*(1-phi_Hs_react[i]),# Col 3  sH -> sh
      (1-e_S_react[i])*(1-g_S_react[i])*phi_Hs_react[i],    # Col 4  sH -> sH
      (1-e_S_react[i])*g_S_react[i]*(1-phi_Hs_react[i]),    # Col 5  sH -> Lh **** Does not match trans matrix--> MCB updated PP 1/18/22 written
      (1-e_S_react[i])*g_S_react[i]*phi_Hs_react[i],        # Col 6  sH -> LH **** Does not match trans matrix--> MCB updated PP 1/18/22 written
      
      # Row 5
      e_L_react[i]*(1-c_H_react[i]),                        # Col 1  Lh -> bh
      e_L_react[i]*c_H_react[i],                            # Col 2  Lh -> bH
      (1-e_L_react[i])*d_L_react[i]*(1-c_H_react[i]),       # Col 3  Lh -> sh **** Does not match trans matrix written --> MCB updated PP 1/18/22
      (1-e_L_react[i])*d_L_react[i]*c_H_react[i],           # Col 4  Lh -> sH **** Does not match trans matrix written --> MCB updated PP 1/18/22
      (1-e_L_react[i])*(1-d_L_react[i])*(1-c_H_react[i]),   # Col 5  Lh -> Lh
      (1-e_L_react[i])*(1-d_L_react[i])*(c_H_react[i]),     # Col 6  Lh -> LH
      
      # Row 6
      e_L_react[i]*(1-phi_HL_react[i]),                      # Col 1  LH -> bh
      e_L_react[i]*phi_HL_react[i],                          # Col 2  LH -> bH
      (1-e_L_react[i])*d_L_react[i]*(1-phi_HL_react[i]),     # Col 3  LH -> sh **** Does not match trans matrix written --> MCB updated PP 1/18/22
      (1-e_L_react[i])*d_L_react[i]*phi_HL_react[i],         # Col 4  LH -> sH **** Does not match trans matrix written --> MCB updated PP 1/18/22
      (1-e_L_react[i])*(1-d_L_react[i])*(1-phi_HL_react[i]), # Col 5  LH -> Lh
      (1-e_L_react[i])*(1-d_L_react[i])*phi_HL_react[i]      # Col 6  LH -> LH
      
    ), 
    nrow = n.states, ncol = n.states, byrow = TRUE)
  }
  
  # Check that all the rows sum to 1
  # rowSums(Bsal_trans_mat_react[, , 1])
  
 
  # At this point, we have set up all the probability/transition matrices 
    # init_conditions
        # Initial conditions
    # trans_mat
        # Transition between occupied and unoccupied host sites pre-Bsal arrival
    # Bsal_trans_mat
        # Transition matrix among 6 states post-Bsal arrival & pre-treatment
    # Bsal_trans_mat_react
        # Transition matrix among 6 states post- Bsal arrival & post-treatment
  # Next, we will simulate the state that each site is in
  
  # As a reminder, this code is testing the:
    # REACTIVE management option
  
  
  # First,
  # Create an empty matrix to store the site state
  z <- array(NA, dim = c(n.sites, n.occasions, n.sim))
  
  
  # Randomly pick the site(s) where Bsal invades
  infected.sites <- sample(1:n.sites, n.Bsal.sites)
  
  
  # Use a for loop to determine what state each site is in
  
  for(n in 1:n.sim){# within each simulation
    
  for(i in 1:n.sites){ # at each site
  
      # First timestep = initial conditions --> #M# determine the status of Bsal/Host at each site at timestep 1
      z[i, 1, n] <- which(rmultinom(1, 1, init_conditions[,n]) == 1)
  }
   
    occ <- length(which(z[,1,n]==2))/n.sites # calculate the proportion of sites occupies at t1
    c_H <- plogis(c_H_alpha + c_H_beta*occ) # function to make host col. dependent on occ. plogis          gets it out of logit scale (inverse logit)
    
    # transition matrix for timestep 1 to 2 which is based on the c_H autologistic function
        trans_mat[,,1,n] <- matrix(c(
          
          # Row 1
          (1-c_H), # Col 1     bh -> bh  - a site remains unoccupied
          c_H,     # Col 2     bh -> bH  - a site is colonized
          
          # Row 2
          (1-phi_Hb[n]), # Col 1     bH -> bh - a site goes extinct
          phi_Hb[n]      # Col 2     bH -> bH - a site remains occupied
          
        ), nrow = 2, ncol = 2, byrow = TRUE)
      
      
      
      # Pre-Bsal invasion timesteps
      for(k in 2:n.pre.yrs){ #M# for each yr between 1st timestep and Bsal
        for(i in 1:n.sites){
        z[i, k, n] <- which(rmultinom(1, 1, trans_mat[z[i, k-1, n], ,k-1, n]) == 1) 
        #M# draw a random number of size 1, based on trans mat probabilities
          # Breaking down the trans_mat indexing:
            # z[i, k-1, n] = the state of site i in time point k-1 in a given simulation
            # then we want all the probabilities in that column (= same timept)
            # And we are on the nth simulation
          # Then, the state at the next time step k goes into z[i,k,n] #M# based on the probabilities in the col?
        }
      
      occ <- length(which(z[,k-1,n]==2))/n.sites # calculate the proportion of sites occupies at the previous timestep for each timestep between first and Bsal arrival
      c_H <- plogis(c_H_alpha + c_H_beta*occ) # function allowing c_H to depend on occ.
        
    #transition matrix for timesteps between 2 and Bsal arrival; based on the c_H autologistic function
        trans_mat[,,k-1,n] <- matrix(c( 
          
          # Row 1
          (1-c_H), # Col 1     bh -> bh  - a site remains unoccupied
          c_H,     # Col 2     bh -> bH  - a site is colonized
          
          # Row 2
          (1-phi_Hb[n]), # Col 1     bH -> bh - a site goes extinct
          phi_Hb[n]      # Col 2     bH -> bH - a site remains occupied
          
        ), nrow = 2, ncol = 2, byrow = TRUE)
        }
      
  }
  
  

  
  ######## Then, Bsal arrives sometime during the last year of n.pre.yrs
    # The dynamics have already occurred and we just switch the state of the site where it arrives
  
  # Determine if the sites that will get Bsal introduced have a host present or absent
    # This will determine if they turn into sh or sH
  
  for(n in 1:n.sim){ # for each simulation
    
    for(i in 1:length(infected.sites)){ # for each infected site
      
      # If a host is absent (bh = state 1) - 
        # the site will convert to sh (state 3)
      if(z[infected.sites[i], n.pre.yrs, n] == 1){
        
        z[infected.sites[i], n.pre.yrs, n] <- 3
        
        }
      
      # If a host is present (bH = state 2) - 
        # the site will convert to sH (state 4)
      if(z[infected.sites[i], n.pre.yrs, n] == 2){
        
        z[infected.sites[i], n.pre.yrs, n] <- 4
        
        }
      
    }
    
  }
  
  # Then, simulate the dynamics - after Bsal arrival
  
  for(n in 1:n.sim){ # within each simulation 
    
    occ_total <- length(which(z[,n.pre.yrs,n]%in% c(2,4,6)))/n.sites # calculate the propor sites occ'd
    c_H <- plogis(c_H_alpha + c_H_beta*occ_total) # function to make host col. depend on occ.
    
    Bsal_trans_mat[,,n.pre.yrs,n] <- matrix(c(
      
      # Row 1
      (1-c_S[n])*(1-c_H),  # Col 1  bh -> bh
      (1-c_S[n])*c_H,      # Col 2  bh -> bH
      c_S[n]*(1-c_H),      # Col 3  bh -> Sh
      c_S[n]*c_H,          # Col 4  bh -> SH
      0,                      # Col 5  bh -> Lh
      0,                      # Col 6  bh -> LH
      
      # Row 2
      (1-c_S[n])*(1-phi_Hb[n]),  # Col 1   bH -> bh #M# when you go H-->h now incorporating persistence likelihoods 
      (1-c_S[n])*phi_Hb[n],      # Col 2   bH -> bH
      c_S[n]*(1-phi_Hb[n]),      # Col 3   bH -> Sh
      c_S[n]*phi_Hb[n],          # Col 4   bH -> SH
      0,                         # Col 5   bH -> Lh
      0,                         # Col 6   bH -> LH
      
      # Row 3
      e_S[n]*(1-c_H),                 # Col 1   Sh -> bh
      e_S[n]*c_H,                     # Col 2   Sh -> bH
      (1-e_S[n])*(1-g_S[n])*(1-c_H),  # Col 3   Sh -> Sh
      (1-e_S[n])*(1-g_S[n])*c_H,      # Col 4   Sh -> SH
      (1-e_S[n])*g_S[n]*(1-c_H),      # Col 5   Sh -> Lh  
      (1-e_S[n])*g_S[n]*c_H,          # Col 6   Sh -> LH  
      
      # Row 4
      e_S[n]*(1-phi_Hs[n]),                 # Col 1  SH -> bh
      e_S[n]*phi_Hs[n] ,                    # Col 2  SH -> bH
      (1-e_S[n])*(1-g_S[n])*(1-phi_Hs[n]),  # Col 3  SH -> Sh
      (1-e_S[n])*(1-g_S[n])*phi_Hs[n],      # Col 4  SH -> SH
      (1-e_S[n])*g_S[n]*(1-phi_Hs[n]),      # Col 5  SH -> Lh 
      (1-e_S[n])*g_S[n]*phi_Hs[n],          # Col 6  SH -> LH 
      
      # Row 5
      e_L[n]*(1-c_H),                 # Col 1   Lh -> bh
      e_L[n]*c_H,                     # Col 2   Lh -> bH
      (1-e_L[n])*d_L[n]*(1-c_H),      # Col 3   Lh -> Sh 
      (1-e_L[n])*(d_L[n])*(c_H),      # Col 4   Lh -> SH 
      (1-e_L[n])*(1-d_L[n])*(1-c_H),  # Col 5   Lh -> Lh
      (1-e_L[n])*(1-d_L[n])*(c_H),    # Col 6   Lh -> LH
      
      # Row 6
      e_L[n]*(1-phi_HL[n]),                 # Col 1  LH -> bh
      e_L[n]*phi_HL[n],                     # Col 2  LH -> bH
      (1-e_L[n])*d_L[n]*(1-phi_HL[n]),      # Col 3  LH -> Sh 
      (1-e_L[n])*d_L[n]*phi_HL[n],          # Col 4  LH -> SH 
      (1-e_L[n])*(1-d_L[n])*(1-phi_HL[n]),  # Col 5  LH -> Lh
      (1-e_L[n])*(1-d_L[n])*phi_HL[n]       # Col 6  LH -> LH
      
    ), 
    nrow = n.states, ncol = n.states, byrow = TRUE)
      
      # for time after Bsal arrival & before reactive treatment
      for(k in (n.pre.yrs+1):n.yrs.post.Bsal.pre.treat){
        
        for(i in 1:n.sites){  # at each site
        z[i, k, n] <- which(rmultinom(1, 1, Bsal_trans_mat[z[i, k-1, n], ,k-1, n]) == 1)
        
      }
      
      occ_total <- length(which(z[,k-1,n]%in% c(2,4,6)))/n.sites # calculate the propor sites occ'd
      c_H <- plogis(c_H_alpha + c_H_beta*occ_total) # function to make host col. depend on occ.
      
      Bsal_trans_mat[,,k-1,n] <- matrix(c(
        
        # Row 1
        (1-c_S[n])*(1-c_H),  # Col 1  bh -> bh, #M# a probability that a site remains without Bsal and without hosts is a product of the prob that a site isn't colonized by Bsal [1-c_S] AND the prob that the site isn't colonized by the host [1-c_H]
        (1-c_S[n])*c_H,      # Col 2  bh -> bH
        c_S[n]*(1-c_H),      # Col 3  bh -> Sh
        c_S[n]*c_H,          # Col 4  bh -> SH
        0,                      # Col 5  bh -> Lh
        0,                      # Col 6  bh -> LH
        
        # Row 2
        (1-c_S[n])*(1-phi_Hb[n]),  # Col 1   bH -> bh #M# when you go H-->h now incorporating persistence likelihoods 
        (1-c_S[n])*phi_Hb[n],      # Col 2   bH -> bH
        c_S[n]*(1-phi_Hb[n]),      # Col 3   bH -> Sh
        c_S[n]*phi_Hb[n],          # Col 4   bH -> SH
        0,                         # Col 5   bH -> Lh
        0,                         # Col 6   bH -> LH
        
        # Row 3
        e_S[n]*(1-c_H),                 # Col 1   Sh -> bh
        e_S[n]*c_H,                     # Col 2   Sh -> bH
        (1-e_S[n])*(1-g_S[n])*(1-c_H),  # Col 3   Sh -> Sh
        (1-e_S[n])*(1-g_S[n])*c_H,      # Col 4   Sh -> SH
        (1-e_S[n])*g_S[n]*(1-c_H),      # Col 5   Sh -> Lh  
        (1-e_S[n])*g_S[n]*c_H,          # Col 6   Sh -> LH  
        
        # Row 4
        e_S[n]*(1-phi_Hs[n]),                 # Col 1  SH -> bh
        e_S[n]*phi_Hs[n] ,                    # Col 2  SH -> bH
        (1-e_S[n])*(1-g_S[n])*(1-phi_Hs[n]),  # Col 3  SH -> Sh
        (1-e_S[n])*(1-g_S[n])*phi_Hs[n],      # Col 4  SH -> SH
        (1-e_S[n])*g_S[n]*(1-phi_Hs[n]),      # Col 5  SH -> Lh 
        (1-e_S[n])*g_S[n]*phi_Hs[n],          # Col 6  SH -> LH 
        
        # Row 5
        e_L[n]*(1-c_H),                 # Col 1   Lh -> bh
        e_L[n]*c_H,                     # Col 2   Lh -> bH
        (1-e_L[n])*d_L[n]*(1-c_H),      # Col 3   Lh -> Sh 
        (1-e_L[n])*(d_L[n])*(c_H),      # Col 4   Lh -> SH 
        (1-e_L[n])*(1-d_L[n])*(1-c_H),  # Col 5   Lh -> Lh
        (1-e_L[n])*(1-d_L[n])*(c_H),    # Col 6   Lh -> LH
        
        # Row 6
        e_L[n]*(1-phi_HL[n]),                 # Col 1  LH -> bh
        e_L[n]*phi_HL[n],                     # Col 2  LH -> bH
        (1-e_L[n])*d_L[n]*(1-phi_HL[n]),      # Col 3  LH -> Sh 
        (1-e_L[n])*d_L[n]*phi_HL[n],          # Col 4  LH -> SH 
        (1-e_L[n])*(1-d_L[n])*(1-phi_HL[n]),  # Col 5  LH -> Lh
        (1-e_L[n])*(1-d_L[n])*phi_HL[n]       # Col 6  LH -> LH
        
      ), 
      nrow = n.states, ncol = n.states, byrow = TRUE)
      }
   }
      
      ##### Treatment applied ######
      
      # Time after Bsal arrival & after reactive treatment  #?# occasions = "timesteps"? noccasion = tottal time steps; n.yrs.post.Bsal.pre.treat = time setsp before treatment
      for(k in (n.yrs.post.Bsal.pre.treat+1):n.occasions){
        
        z[i, k, n] <- which(rmultinom(1, 1, Bsal_trans_mat_react[z[i, k-1, n], , n]) ==1)
        
      }
      
    }
    
  
  
  
  # By this point, we have the site states, and we need to summarize them
  
  # Create an empty matrix to summarize the proportion of sites in each state at each time step
  Bsal_sum <- array(NA, dim = c(n.states, ncol(z), n))
  
  # Add rownames
  rownames(Bsal_sum) <- c("bh", "bH", "sh", "sH", "Lh", "LH")
  
  # Loop through each timestep and state to determine the proportion of sites in each state
  for(n in 1:n.sim){ # for each sim
    
    for(i in 1:n.occasions){ # at each timestep
      
      for(k in 1:n.states){ # for each state
        
        # Count up the number of sites in each state and divide by the total number of sites to obtain a proportion
        Bsal_sum[k, i, n] <- length(which(z[, i, n] == k)) / n.sites
        
      }
      
    }
    
  }
  

  # Convert the matrix into long format
  Bsal_sum.long <- melt(Bsal_sum)
  
  # Add column names
  colnames(Bsal_sum.long) <- c("state", "n.occasion", "n.sim" ,"proportion")
  
  
  # Create a new object with just the data you want
    # We want the following information for the last occasion 
      # state
      # proportion of sites
      # simulation ID
  #M# so giving the final "status" of each site after 10 yrs across all the sims #?#
  last.occ <- data.frame(state = Bsal_sum.long$state[which(Bsal_sum.long$n.occasion == n.occasions)],
                         proportions = Bsal_sum.long$proportion[which(Bsal_sum.long$n.occasion == n.occasions)],
                         n.sim = Bsal_sum.long$n.sim[which(Bsal_sum.long$n.occasion == n.occasions)]) #?# which statements? pulls out the statuses at the last occastion/timesteps
  
  # Then add a new column and summarize the data 
    # Summarize the mean + 95% CI proportion of hosts in each state on the last occasion
  params_dat <- last.occ %>%
    add_column(host = rep(c("h", "H", "h", "H", "h", "H"), times = n.sim))%>% #?# This is just replicating the host pres/abs for each sim which is now in long format? -- Yes
    group_by(state, host) %>%
    summarise(mean = mean(proportions),
              lower = quantile(proportions, prob = c(0.025, 0.975))[1],
              upper = quantile(proportions, prob = c(0.025, 0.975))[2])
  
  # Objects to return
  return(list(params_dat = params_dat,      # Summarized data
              trajectory = Bsal_sum.long))  # Raw data from simulation
  
}





# 5. Run the different scenarios ---------------------------------------------------------------



# Total number of timesteps 
tot.time <- 40


# 2 timesteps
pre.val <- 2 + n.pre.yrs  #?# this means reactive management only after 2 years? -- yes
post.val <- tot.time - pre.val
params_dat_2postBsalarrival <- Bsal.sim.react(n.yrs.post.Bsal.pre.treat = pre.val,
                                              n.post.Bsal.post.treat = post.val,
                                              n.occasions = n.pre.yrs + pre.val + post.val)


# 4 timesteps
pre.val <- 4 + n.pre.yrs
post.val <- tot.time - pre.val
params_dat_4postBsalarrival <- Bsal.sim.react(n.yrs.post.Bsal.pre.treat = pre.val,
                                              n.post.Bsal.post.treat = post.val,
                                              n.occasions = n.pre.yrs + pre.val + post.val)

# 8 timesteps
pre.val <- 8 + n.pre.yrs
post.val <- tot.time - pre.val
params_dat_8postBsalarrival <- Bsal.sim.react(n.yrs.post.Bsal.pre.treat = pre.val,
                                              n.post.Bsal.post.treat = post.val,
                                              n.occasions = n.pre.yrs + pre.val + post.val)


# 16 timesteps
pre.val <- 16 + n.pre.yrs
post.val <- tot.time - pre.val
params_dat_16postBsalarrival <- Bsal.sim.react(n.yrs.post.Bsal.pre.treat = pre.val,
                                               n.post.Bsal.post.treat = post.val,
                                               n.occasions = n.pre.yrs + pre.val + post.val)




# 6. Create a figure with % of hosts in each group over time by scenario (#?# = when reactive treatment is done?) ---------------------------------------------------------------


# Create one dataframe with all the simulations
vals <- rbind(as.data.frame(cbind(yrs = 2, params_dat_2postBsalarrival$trajectory)),
              as.data.frame(cbind(yrs = 4, params_dat_4postBsalarrival$trajectory)),
              as.data.frame(cbind(yrs = 8, params_dat_8postBsalarrival$trajectory)),
              as.data.frame(cbind(yrs = 16, params_dat_16postBsalarrival$trajectory))
)

# Change column names
colnames(vals) <- c("react_time", "state", "n.occasion", "n.sim" , "proportion")

# Create dataframe with lines - will be used for the ggplot figure
react_time <- data.frame(state = rep(c("bh", "bH", "sh", "sH", "Lh", "LH"), times = 4),
                         react_time = rep(c(2, 4, 8, 16), each = 6),
                         lines = rep(c(2, 4, 8, 16) + n.pre.yrs, each = 6))


# Make the figure
ggplot(data = vals, aes(x = n.occasion, y = proportion)) + 
  geom_line(aes(col = as.factor(n.sim),alpha = 0.5))+
  facet_grid(state ~ react_time)+
  xlab("Time step")+
  ylab("Proportion of sites")+
  theme_bw()+
  geom_vline(xintercept = n.pre.yrs, lty = 2)+
  geom_vline(data = react_time, aes(xintercept = lines), lty = 2, col = "red")+
  geom_smooth()+
  theme(axis.title = element_text(size = 17),
        axis.text =  element_text(size = 17),
        strip.text =  element_text(size = 17),
        legend.position = "none")

# Save the figure
# ggsave("./Figures/React_sims_proportion_vs_time_2020_11_16.pdf", height = 12, width = 13)





# 7. Create a summary of % of hosts remaining in each group ---------------------------------------------------------------


# Create one dataframe
vals <- rbind(as.data.frame(cbind(yrs = 2, params_dat_2postBsalarrival$params_dat)),
              as.data.frame(cbind(yrs = 4, params_dat_4postBsalarrival$params_dat)),
              as.data.frame(cbind(yrs = 8, params_dat_8postBsalarrival$params_dat)),
              as.data.frame(cbind(yrs = 16, params_dat_16postBsalarrival$params_dat))
)

# Change column names
colnames(vals) <- c("yrs", "State", "host","mean", "lower", "upper")

# Add a column with host persisting or extinct
vals$Host <- ifelse(vals$host == "H", "Host persisting", "Host extinct")

# Make the figure
ggplot(data = vals, aes(x = yrs, y = mean, col = State)) + 
  geom_line(position = position_dodge(width = 0.9))+
  geom_pointrange(aes(x = yrs, ymin = lower, ymax = upper), position = position_dodge(width = 0.9))+
  facet_wrap(~Host)+
  xlab("Number of timesteps after Bsal arrival and before treatment")+
  ylab("Proportion of sites")+
  theme_bw()+
  theme(axis.title = element_text(size = 17),
        axis.text =  element_text(size = 17),
        strip.text =  element_text(size = 17))

##?## constraint that makes more biologically realistic; bounce back is quickly
# deterministic model reaches equilibrium quickly; stochastic - more variability 

# Save the figure
# ggsave("./Figures/React_sims_yrs_state.pdf", height = 6, width = 8)


# End script

