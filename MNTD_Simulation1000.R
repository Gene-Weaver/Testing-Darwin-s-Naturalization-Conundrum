# Code for: https://dx.doi.org/10.1111/ddi.12861
#
#     Ng, J., W. Weaver,* and R.G. Laport. In Press. 
#         Testing Darwin's naturalization conundrum using phylogenetic relationships: generalizable 
#         patterns across disparate communities. Diversity and Distributions. 

###### Set the directory location for output files. If it is in the working directory use ""
DIRECTORY <- ""

# ALL_SITES_OUTPUT columns 2-6 ==> MNTD values for TREE_Obs
#                  columns 7-21 ==> P-values where null distribution is from the observed MNTDs for the 1000resampledtrees
#                                   Can be used to see how TREE_Obs differs from the 1000resampledtrees
#                  columns 22-36 ==> P-values where null distribution is from the randomized invasive 1000resampledtrees
#                                   Can be used as a statistical null dist. for TREE_Obs

library(ape)
library(picante)
library(dplyr)
library(tictoc)
library(evobiR)
library(phyloch)

source("./Will/0_Sandbox/MNTD_MOD.R")
source("./Will/0_Sandbox/SES_MNTD_MOD.R")
source("./Will/9_Ecology/Invasive_Binomial_Builder.R")
source("./Will/0_Sandbox/Status.0.1.Builder.R")

### Trees
TREE_Obs <- read.nexus("./NEON_MCC_wgymnoclade.tre")
TREE_1000 <- read.nexus("./NEONtree_wgymnoclade_1000.trees")

### Abundance Data
ABUND_RAW <- read.csv("./Abundance_Data_Species_as_Col_SEVEN.csv") #326 species
ABUND_ALL_SITES <- ABUND_RAW[,which(colnames(ABUND_RAW)%in%TREE_Obs$tip.label)] #319 species

### Get Site Names
SITENAMES <- data.frame(as.character(ABUND_RAW[,1]))
colnames(SITENAMES) <- "Site"

### Import USDA Invasive csv
INVASIVE_0 <- read.csv("./USDA_Invasive.csv")
### Invasive Binomial Builder
temp <- Invasive.Binomial.Builder(INVASIVE_0)
IBinomial <- data.frame(temp$IBinomial)
INVASIVE <- data.frame(temp$Invasive)
colnames(IBinomial) <- c("IBinomial")
colnames(INVASIVE) <- c("Invasive")

### List of all invasives at all NEON sites
NEON_INVASIVE_SPECIES <- ABUND_RAW[,which(colnames(ABUND_RAW)%in%IBinomial$IBinomial)]

### Regional Species Pools
RSP_HARV <- read.csv("./Regional_Species_Pool_WP_HARV.csv")
RSP_JERC <- read.csv("./Regional_Species_Pool_WP_JERC.csv")
RSP_ORNL <- read.csv("./Regional_Species_Pool_WP_ORNL.csv")
RSP_OSBS <- read.csv("./Regional_Species_Pool_WP_OSBS.csv")
RSP_SCBI <- read.csv("./Regional_Species_Pool_WP_SCBI.csv")
RSP_SERC <- read.csv("./Regional_Species_Pool_WP_SERC.csv")
RSP_TALL <- read.csv("./Regional_Species_Pool_WP_TALL.csv")
RSP_ALL_SITES <- list(RSP_HARV,RSP_JERC,RSP_ORNL,RSP_OSBS,RSP_SCBI,RSP_SERC,RSP_TALL)
RSP_NEON <- read.csv("./Regional_Species_Pool_WP_NEON.csv")

print(setdiff(IGNORE[,1],colnames(ABUND_ALL_SITES)))
AAA <- data.frame(colnames(ABUND_ALL_SITES))
### List of species that should not appear in randomizations + one species that was duplicated
IGNORE <- data.frame(c("Macfadyena_unguis.cati","Smilax_bona.nox","Abies_alba",
                       "Pinus_thunbergii","Callitris_glaucophylla","Platycladus_orientalis",
                       "Taxus_baccata","Taxus_cuspidata","Lygodium_microphyllum",
                       "Lygodium_japonicum"))
### Drop tips from Obs tree, drop tips for 1000 tree has to occur when it's used
TREE_Obs <- drop.tip2(TREE_Obs,as.character(IGNORE[,1]))
### MNTD Function
#Calc.MNTD(ABUND_ALL_SITES,IBinomial$IBinomial,TREE = "TREE-HERE",n_runs = 999)
#ABUND <- ABUND_SIM
#InvasiveBinomial <- IBinomial$IBinomial
#TREE <- SIMTREE
#n_runs <- 999


#MNTD_OBSERVED_wRSP <- Calc.MNTD(RSP_AND_SITE_formatted,IBinomial$IBinomial,TREE = TREE_SITE_wRSP,n_runs = 999)
# ABUND <- RSP_AND_SITE_formatted
# InvasiveBinomial <- IBinomial$IBinomial
# TREE <- TREE_SITE_wRSP
# # <- Calc.MNTD(ABUND_SITE,IBinomial$IBinomial,TREE = TREE_SITE,n_runs = 999)
# ABUND <- ABUND_SITE
# InvasiveBinomial <- IBinomial$IBinomial
# TREE <- TREE_SITE


Calc.MNTD <- function(ABUND,InvasiveBinomial,TREE,n_runs){
    #Build Components 
    temp_out <- Status.0.1.Builder(ABUND,InvasiveBinomial)
    COMM <- temp_out[[1]]
    COMM_ALL <- rbind(COMM,matrix(data=1,nrow = 1,ncol=(length(COMM)/2)))
    #SES.MNTD
    MNTD_OUT <- ses.mntd(COMM_ALL, cophenetic(TREE), abundance.weighted=FALSE,runs = n_runs)
    #Parse SES.MNTD
    MNTD_ALL_species_at_site <- MNTD_OUT[3,2]
    MNTD_between_NATIVE_species_at_site <- MNTD_OUT[2,2]
    MNTD_between_NONNATIVE_species_at_site <- MNTD_OUT[1,2]
    #SES.MNTD.MOD
    MNTD_MOD_OUT <- SES_MNTD_MOD(COMM_ALL, cophenetic(TREE), abundance.weighted=FALSE,runs = n_runs)
    #Parse SES.MNTD.MOD
    #MNTD_MOD_ALL_species_at_site <- MNTD_MOD_OUT[3,2]
    MNTD_MOD_NN_N <- MNTD_MOD_OUT[2,2] #Should be smaller (most times)
    MNTD_MOD_N_NN <- MNTD_MOD_OUT[1,2] #Should be a bigger number
    return(cbind.data.frame(MNTD_ALL_species_at_site,MNTD_between_NATIVE_species_at_site,MNTD_between_NONNATIVE_species_at_site,
                            MNTD_MOD_NN_N,MNTD_MOD_N_NN))
    #return(cbind.data.frame(MNTD_ALL_species_at_site,MNTD_between_NATIVE_species_at_site))
}

### Sim Function
Build.Simulation <- function(tree,SITE_INVASIVE_COUNT,RSP_INVASIVE_SPECIES,tree_names_n_only){
    #BUILD SIM TREE
    SIM_TREE <- tree
    DROP <- setdiff(colnames(RSP_INVASIVE_SPECIES),SIM_TREE$tip.label)
    SIM_TREE <- drop.tip(SIM_TREE,as.character(DROP))
    #print(length(colnames(RSP_INVASIVE_SPECIES)))
    RSP_INVASIVE_SPECIES <- RSP_INVASIVE_SPECIES[,which(colnames(RSP_INVASIVE_SPECIES)%in%IGNORE[,1]==FALSE)]
    RSP_INVASIVE_SPECIES <- RSP_INVASIVE_SPECIES[,which(colnames(RSP_INVASIVE_SPECIES)%in%SIM_TREE$tip.label==TRUE)]
    #print(length(colnames(RSP_INVASIVE_SPECIES)))
    #print("...")
    # take random sample of invasive species and add the names to sim_species 
    sim_I_SAMPLE_0 <- data.frame(sample(colnames(RSP_INVASIVE_SPECIES),SITE_INVASIVE_COUNT,replace = FALSE))
    
    # sim_I_SAMPLE_00 <- sim_I_SAMPLE_0
    # colnames(sim_I_SAMPLE_00) <- "I_Species"
    # print(sim_I_SAMPLE_00)
    
    #Native names carried through simulations per site
    sim_N <- data.frame(tree_names_n_only)
    
    #blank dataframe to fill with random invasive species
    sim_I_SAMPLE <- rep(NA,SITE_INVASIVE_COUNT)
    
    #to get correct formatting
    for(i in 1:SITE_INVASIVE_COUNT){
        sim_I_SAMPLE[i] <- as.character(sim_I_SAMPLE_0[i,1])
    }
    sim_I_SAMPLE <- data.frame(sim_I_SAMPLE)
    #print(sim_I_SAMPLE)
    
    colnames(sim_I_SAMPLE) <- "Sim_Species_All"
    sim_N <- data.frame(as.character(sim_N[1,]))
    colnames(sim_N) <- "Sim_Species_All"
    
    #merge sim_I and natives
    sim_species <- rbind.data.frame(sim_I_SAMPLE[1],sim_N[1])
    
    #print(sim_species)
    
    #### drop all other tips to make it site specific with the random Invasives 
    sim_tree_RAND <- drop.tip(SIM_TREE,which(SIM_TREE$tip.label%in%sim_species$Sim_Species_All==FALSE))
    
    #print(length(sim_tree_RAND$tip.label))
    
    return(sim_tree_RAND)
}

Calc.P.Value <- function(OBSERVED,SIMULATION_MNTD,numReps){
    p_v <- sum( SIMULATION_MNTD > OBSERVED )/numReps
    # print("P.Value")
    # print(paste(p_v,SIMULATION_MNTD,OBSERVED,numReps))
    if((p_v >0.5) == TRUE){
      p_v <- (1-p_v)
    }
    d <- density(SIMULATION_MNTD)
    q2_5 <- quantile(SIMULATION_MNTD, 0.025)
    q97_5 <- quantile(SIMULATION_MNTD, 0.975)
    #x1 <- min(d$x[which(d$x > q2_5)])
    #x2 <- max(d$x[which(d$x < q97_5)])
    return(p_v)
}

Calc.P.Value.Details <- function(OBSERVED,SIMULATION_MNTD,numReps){
    p_v <- sum( SIMULATION_MNTD > OBSERVED )/numReps
    # print("P.Value.Detail")
    # print(paste(p_v,SIMULATION_MNTD,OBSERVED,numReps))
    if((p_v >0.5) == TRUE){
      p_v <- (1-p_v)
    }
    d <- density(SIMULATION_MNTD)
    q2_5 <- quantile(SIMULATION_MNTD, 0.025)
    q97_5 <- quantile(SIMULATION_MNTD, 0.975)
    #x1 <- min(d$x[which(d$x > q2_5)])
    #x2 <- max(d$x[which(d$x < q97_5)])
    return(list(p_v,q2_5,q97_5))
}

# i = 1
# j = 1
# k = 1
n <- 1000
ALL_SITES_OUTPUT <- data.frame()
for(i in 1:length(SITENAMES[,1])){
    tic("Site")
    #SITE_OUTPUT <- data.frame()
    CURRENT_SITE <- SITENAMES[i,]
    print(CURRENT_SITE)
    
    ### RSP
    RSP_SITE <- data.frame(RSP_ALL_SITES[[i]])
    colnames(RSP_SITE) <- "Species"
    RSP_SITE_N <- data.frame(RSP_SITE[which(RSP_SITE[,1]%in%IBinomial$IBinomial==FALSE),1])
    RSP_SITE_NN <- data.frame(RSP_SITE[which(RSP_SITE[,1]%in%IBinomial$IBinomial==TRUE),1])
    RSP_SITE_NN_formatted <- t(RSP_SITE_NN)
    colnames(RSP_SITE_NN_formatted) <- RSP_SITE_NN[,1]
    row.names(RSP_SITE_NN_formatted) <- NULL
    RSP_SITE_NN_formatted <- data.frame(RSP_SITE_NN_formatted)
    
    
    ### Prune Tree to Site
    # *** SPECIES_SITE is used for no RSP, RSP_AND_SITE_formatted for RSP
    
    ABUND_SITE_ALL <- rbind(colnames(ABUND_ALL_SITES),ABUND_ALL_SITES[CURRENT_SITE,])
    SPECIES_SITE <- ABUND_SITE_ALL[1,][which(ABUND_SITE_ALL[2,]>=1)]
    SPECIES_SITE_formatted <- data.frame(t(SPECIES_SITE))
    colnames(SPECIES_SITE_formatted)  <- "Species"
    
    RSP_AND_SITE <- rbind(SPECIES_SITE_formatted,RSP_SITE)
    RSP_AND_SITE <- data.frame(unique(RSP_AND_SITE[,1]))
    row.names(RSP_AND_SITE) <- NULL
    RSP_AND_SITE_formatted <- data.frame(t(RSP_AND_SITE))
    colnames(RSP_AND_SITE_formatted) <- RSP_AND_SITE[,1]
    row.names(RSP_AND_SITE_formatted) <- NULL
    
    # Add pseudo abundance data, all = 1
    add <- data.frame(t(rep("1",length(RSP_AND_SITE[,1]))))
    colnames(add) <- RSP_AND_SITE[,1]
    row.names(add) <- NULL
    RSP_AND_SITE_formatted <- rbind.data.frame(RSP_AND_SITE_formatted,add)
    
    # Prune Trees
    TREE_SITE_wRSP <- drop.tip2(TREE_Obs,which(TREE_Obs$tip.label%in%colnames(RSP_AND_SITE_formatted)==FALSE))
    TREE_SITE <- drop.tip2(TREE_Obs,which(TREE_Obs$tip.label%in%SPECIES_SITE==FALSE))
    
    # Remove RSP Species that are not in the tree
    RSP_AND_SITE_formatted <- RSP_AND_SITE_formatted[,which(colnames(RSP_AND_SITE_formatted)%in%TREE_SITE_wRSP$tip.label==TRUE)]
    
    ### Prune Abundance to Site
    ABUND_SITE <- ABUND_SITE_ALL[,which(colnames(ABUND_SITE_ALL)%in%TREE_SITE$tip.label)]
    
    ### Stats for Summary Tree
    MNTD_OBSERVED <- Calc.MNTD(ABUND_SITE,IBinomial$IBinomial,TREE = TREE_SITE,n_runs = 999)
    
    MNTD_OBSERVED_wRSP <- Calc.MNTD(RSP_AND_SITE_formatted,IBinomial$IBinomial,TREE = TREE_SITE_wRSP,n_runs = 999)
    
    # MNTD_OBSERVED_Between_ALL <- MNTD_OBSERVED[1]
    # MNTD_OBSERVED_Between_NATIVE <- MNTD_OBSERVED[2]
    # MNTD_OBSERVED_Between_NONNATIVE <- MNTD_OBSERVED[3]
    # MNTD_MOD_OBSERVED_Between_NN-N <- MNTD_OBSERVED[4]
    # MNTD_MOD_OBSERVED_Between_N_NN <- MNTD_OBSERVED[5]
    
    ### n Invasive
    I_COUNT <- length(which(colnames(SPECIES_SITE)%in%IBinomial$IBinomial==TRUE))
    tree_names_n_only <- ABUND_SITE[1,which(colnames(ABUND_SITE)%in%IBinomial$IBinomial==FALSE)]
    
    SUBTREE_MNTD_OBSERVED_Output <- data.frame()
    SUBTREE_MNTD_RAND_Output <- data.frame()
    for(j in 1:length(TREE_1000)){
        tic("Subtree")
        print(paste(CURRENT_SITE," ",j," / ",length(TREE_1000)))
        SUBTREE <- TREE_1000[[j]]
        
        SUBTREE_SITE <- drop.tip2(SUBTREE,which(SUBTREE$tip.label%in%SPECIES_SITE==FALSE))
        SUBTREE_SITE <- drop.tip(SUBTREE_SITE,as.character(IGNORE[,1]))
        SUBTREE_SITE_wRSP <- drop.tip2(SUBTREE,which(SUBTREE$tip.label%in%colnames(RSP_AND_SITE_formatted)==FALSE))
        SUBTREE_SITE_wRSP <- drop.tip(SUBTREE_SITE_wRSP,as.character(IGNORE[,1]))
        
        ### 1000 Randomizations for SUBTREE ***OUTPUT => one line wi th quartiles, means, medians, p-value vs. SUBTREE_Obs
        SUBTREE_Rand_Output <- data.frame(NA,NA,NA,NA,NA)
        colnames(SUBTREE_Rand_Output) <- c("MNTD_ALL_species_at_site","MNTD_between_NATIVE_species_at_site","MNTD_between_NONNATIVE_species_at_site","MNTD_MOD_NN_N","MNTD_MOD_N_NN")
        for(k in 1:n){
            SIMTREE <- Build.Simulation(SUBTREE_SITE_wRSP,I_COUNT,RSP_SITE_NN_formatted,tree_names_n_only)
            if(length(SIMTREE$tip.label) != length(SUBTREE_SITE$tip.label)){
              print(paste(length(SUBTREE_SITE$tip.label)," <- Obs...Sim -> ",length(SIMTREE$tip.label),sep=""))
            }
            
            ABUND_SIM <- RSP_AND_SITE_formatted[,colnames(RSP_AND_SITE_formatted)%in%SIMTREE$tip.label==TRUE]
            SIMTREE_MNTD <- data.frame(Calc.MNTD(ABUND_SIM,IBinomial$IBinomial,TREE = SIMTREE,n_runs = 1))
            colnames(SIMTREE_MNTD) <- c("MNTD_ALL_species_at_site",
                                        "MNTD_between_NATIVE_species_at_site",
                                        "MNTD_between_NONNATIVE_species_at_site",
                                        "MNTD_MOD_NN_N",
                                        "MNTD_MOD_N_NN")
            SUBTREE_Rand_Output <- rbind.data.frame(SUBTREE_Rand_Output,SIMTREE_MNTD)
        }#END Randomization Loop
        SUBTREE_Rand_Output <- SUBTREE_Rand_Output[-1,]   ###Don't need to keep these values, just use to get p-values
        
        ### Stats for SUBTREE
        SUBTREE_MNTD_OBSERVED <- data.frame(Calc.MNTD(ABUND_SITE,IBinomial$IBinomial,TREE = SUBTREE_SITE,n_runs = 999))
        
        ### Randomizations for 1000subtrees to generate null for TREE_Obs
        SUBTREE_RAND <- Build.Simulation(SUBTREE,I_COUNT,RSP_SITE_NN_formatted,tree_names_n_only)
        ABUND_SIM <- RSP_AND_SITE_formatted[,colnames(RSP_AND_SITE_formatted)%in%SUBTREE_RAND$tip.label==TRUE]
        SUBTREE_MNTD_RAND <- data.frame(Calc.MNTD(ABUND_SIM,IBinomial$IBinomial,TREE = SUBTREE_RAND,n_runs = 999))
        
        # SUBTREE_MNTD_OBSERVED_Between_ALL <- SUBTREE_MNTD_OBSERVED[1]
        # SUBTREE_MNTD_OBSERVED_Between_NATIVE <- SUBTREE_MNTD_OBSERVED[2] etc...
        
        ### P-Values
        P_ALL <- Calc.P.Value(SUBTREE_MNTD_OBSERVED[,1],SUBTREE_Rand_Output[,1],n)
        P_NATIVE <- Calc.P.Value(SUBTREE_MNTD_OBSERVED[,2],SUBTREE_Rand_Output[,2],n)
        P_NONNATIVE <- Calc.P.Value(SUBTREE_MNTD_OBSERVED[,3],SUBTREE_Rand_Output[,3],n)
        P_NN_N <- Calc.P.Value(SUBTREE_MNTD_OBSERVED[,4],SUBTREE_Rand_Output[,4],n)
        P_N_NN <-Calc.P.Value(SUBTREE_MNTD_OBSERVED[,5],SUBTREE_Rand_Output[,5],n)
        
        SUBTREE_MNTD_RAND_Output <- rbind.data.frame(SUBTREE_MNTD_RAND_Output,SUBTREE_MNTD_RAND)
        SUBTREE_MNTD_OBSERVED_Output <- rbind.data.frame(SUBTREE_MNTD_OBSERVED_Output,
                                                         cbind.data.frame(j,SUBTREE_MNTD_OBSERVED,
                                                                          P_ALL,P_NATIVE,P_NONNATIVE,P_NN_N,P_N_NN))
        toc()
    }#END TREE_1000 Loop
    
    ######################## For observed values in 1000resampledtrees
    colnames(SUBTREE_MNTD_OBSERVED_Output) <- c("SUBTREE","MNTD_ALL_species_at_site","MNTD_between_NATIVE_species_at_site",
                                                "MNTD_between_NONNATIVE_species_at_site","MNTD_MOD_NN_N","MNTD_MOD_N_NN",
                                                "P_ALL","P_NATIVE","P_NONNATIVE","P_NN_N","P_N_NN")
    filename_SUBTREE_MNTD_OBSERVED_Output <- paste("./",DIRECTORY,CURRENT_SITE,"_MNTD_Observed1000Trees_7-5-18.csv",sep="")
    write.csv(SUBTREE_MNTD_OBSERVED_Output,filename_SUBTREE_MNTD_OBSERVED_Output,row.names=FALSE)
    ### 1000 Tree p-value for using observed mntd(1000trees)
    P_ALL <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,1],SUBTREE_MNTD_OBSERVED_Output[,2],n))
    P_NATIVE <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,2],SUBTREE_MNTD_OBSERVED_Output[,3],n))
    P_NONNATIVE <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,3],SUBTREE_MNTD_OBSERVED_Output[,4],n))
    P_NN_N <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,4],SUBTREE_MNTD_OBSERVED_Output[,5],n))
    P_N_NN <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,5],SUBTREE_MNTD_OBSERVED_Output[,6],n))
    colnames(P_ALL) <- c("P_Value-All","q2_5-All","q97_5-All")
    colnames(P_NATIVE) <- c("P_Value-Native","q2_5-Native","q97_5-Native")
    colnames(P_NONNATIVE) <- c("P_Value-NonNative","q2_5-NonNative","q97_5-NonNative")
    colnames(P_NN_N) <- c("P_Value-NN_N","q2_5-NN_N","q97_5-NN_N")
    colnames(P_N_NN) <- c("P_Value-N_NN","q2_5-N_NN","q97_5-N_NN")
    
    ######################## For randomized invasives in 1000resampledtrees
    colnames(SUBTREE_MNTD_RAND_Output) <- c("MNTD_ALL_species_at_site_RAND",
                                            "MNTD_between_NATIVE_species_at_site_RAND",
                                            "MNTD_between_NONNATIVE_species_at_site_RAND",
                                            "MNTD_NN_N_species_at_site_RAND",
                                            "MNTD_N_NN_species_at_site_RAND")
    filename_SUBTREE_MNTD_RAND_Output <- paste("./",DIRECTORY,CURRENT_SITE,"_MNTD_RandomizedInvasives1000Trees_7-5-18.csv",sep="")
    write.csv(SUBTREE_MNTD_RAND_Output,filename_SUBTREE_MNTD_RAND_Output,row.names=FALSE)
    ### 1000 Tree p-value for using RANDOMIZED INVASIVES mntd(1000trees)
    P_ALL_RI <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,1],SUBTREE_MNTD_RAND_Output[,1],n))
    P_NATIVE_RI <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,2],SUBTREE_MNTD_RAND_Output[,2],n))
    P_NONNATIVE_RI <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,3],SUBTREE_MNTD_RAND_Output[,3],n))
    P_NN_N_RI <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,4],SUBTREE_MNTD_RAND_Output[,4],n))
    P_N_NN_RI <- data.frame(Calc.P.Value.Details(MNTD_OBSERVED[,5],SUBTREE_MNTD_RAND_Output[,5],n))
    colnames(P_ALL_RI) <- c("P_Value-All_RI","q2_5-All_RI","q97_5-All_RI")
    colnames(P_NATIVE_RI) <- c("P_Value-Native_RI","q2_5-Native_RI","q97_5-Native_RI")
    colnames(P_NONNATIVE_RI) <- c("P_Value-NonNative_RI","q2_5-NonNative_RI","q97_5-NonNative_RI")
    colnames(P_NN_N_RI) <- c("P_Value-NN_N_RI","q2_5-NN_N_RI","q97_5-NN_N_RI")
    colnames(P_N_NN_RI) <- c("P_Value-N_NN_RI","q2_5-N_NN_RI","q97_5-N_NN_RI")
    
    ### Add to ALL_SITES_OUTPUT
    # MNTD_OBSERVED_Between_ALL <- MNTD_OBSERVED[1]
    # MNTD_OBSERVED_Between_NATIVE <- MNTD_OBSERVED[2] etc...
    ALL_SITES_OUTPUT <- rbind.data.frame(ALL_SITES_OUTPUT,cbind.data.frame(CURRENT_SITE,MNTD_OBSERVED[1],MNTD_OBSERVED[2],MNTD_OBSERVED[3],MNTD_OBSERVED[4],MNTD_OBSERVED[5],
                                                                           P_ALL,P_NATIVE,P_NONNATIVE,P_NN_N,P_N_NN,
                                                                           P_ALL_RI,P_NATIVE_RI,P_NONNATIVE_RI,P_NN_N_RI,P_N_NN_RI,
                                                                           MNTD_OBSERVED_wRSP[1],MNTD_OBSERVED_wRSP[2],MNTD_OBSERVED_wRSP[3],MNTD_OBSERVED_wRSP[4],MNTD_OBSERVED_wRSP[5]))
    toc()
}#END Site Loop

### Write ALL_SITES_OUTPUT
filename_ALL_SITES_OUTPUT <- paste("./",DIRECTORY,"NEON_MNTD_Observed.csv",sep="")
write.csv(ALL_SITES_OUTPUT,filename_ALL_SITES_OUTPUT,row.names=FALSE)



