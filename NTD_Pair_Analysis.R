# Code for: https://dx.doi.org/10.1111/ddi.12861
#
#     Ng, J., W. Weaver,* and R.G. Laport. In Press. 
#         Testing Darwin's naturalization conundrum using phylogenetic relationships: generalizable 
#         patterns across disparate communities. Diversity and Distributions. 

###### Set the directory location for output files. If it is in the working directory use ""
DIRECTORY <- ""

# NTD Pair Analysis


library(ape)
library(picante)
library(dplyr)
library(tictoc)
require(plyr)
source("./Toolbox/MNTD_MOD.R")
source("./Toolbox/SES_MNTD_MOD.R")
source("./Toolbox/Invasive_Binomial_Builder.R")
source("./Toolbox/Status.0.1.Builder.R")

### Trees
TREE_Obs <- read.nexus("./NEON_MCC_wgymnoclade_NEW.tre")

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

### Regional Species Pools
RSP_HARV <- read.csv("./Will/Revisions/Regional_Species_Pool_WP_HARV.csv")
RSP_JERC <- read.csv("./Will/Revisions/Regional_Species_Pool_WP_JERC.csv")
RSP_ORNL <- read.csv("./Will/Revisions/Regional_Species_Pool_WP_ORNL.csv")
RSP_OSBS <- read.csv("./Will/Revisions/Regional_Species_Pool_WP_OSBS.csv")
RSP_SCBI <- read.csv("./Will/Revisions/Regional_Species_Pool_WP_SCBI.csv")
RSP_SERC <- read.csv("./Will/Revisions/Regional_Species_Pool_WP_SERC.csv")
RSP_TALL <- read.csv("./Will/Revisions/Regional_Species_Pool_WP_TALL.csv")
RSP_ALL_SITES <- list(RSP_HARV,RSP_JERC,RSP_ORNL,RSP_OSBS,RSP_SCBI,RSP_SERC,RSP_TALL)
RSP_NEON <- read.csv("./Will/Revisions/Regional_Species_Pool_WP_NEON.csv")

IGNORE <- data.frame(c("Macfadyena_unguis.cati","Smilax_bona.nox","Abies_alba",
                       "Pinus_thunbergii","Callitris_glaucophylla","Platycladus_orientalis",
                       "Taxus_baccata","Taxus_cuspidata","Lygodium_microphyllum",
                       "Lygodium_japonicum"))
TREE_Obs <- drop.tip(TREE_Obs,as.character(IGNORE[,1]))

### List of all invasives at all NEON sites
NEON_INVASIVE_SPECIES <- ABUND_RAW[,which(colnames(ABUND_RAW)%in%IBinomial$IBinomial)]

### MNTD Function
#Calc.MNTD(ABUND_ALL_SITES,IBinomial$IBinomial,TREE = "TREE-HERE",n_runs = 999)
#ABUND <- ABUND_SIM
#InvasiveBinomial <- IBinomial$IBinomial
#TREE <- SIMTREE
#n_runs <- 999
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
    if(p_v >.5){p_v <- 1-p_v}
    d <- density(SIMULATION_MNTD)
    q2_5 <- quantile(SIMULATION_MNTD, 0.05)
    q97_5 <- quantile(SIMULATION_MNTD, 0.95)
    #x1 <- min(d$x[which(d$x > q2_5)])
    #x2 <- max(d$x[which(d$x < q97_5)])
    return(p_v)
}

Calc.P.Value.Details <- function(OBSERVED,SIMULATION_MNTD,numReps){
    p_v <- sum( SIMULATION_MNTD > OBSERVED )/numReps
    if(p_v >.5){p_v <- 1-p_v}
    d <- density(SIMULATION_MNTD)
    q2_5 <- quantile(SIMULATION_MNTD, 0.05)
    q97_5 <- quantile(SIMULATION_MNTD, 0.95)
    #x1 <- min(d$x[which(d$x > q2_5)])
    #x2 <- max(d$x[which(d$x < q97_5)])
    return(list(p_v,q2_5,q97_5))
}

MNTD_MOD_return_species_names <- function (samp, dis, abundance.weighted = FALSE){
    N <- dim(samp)[1]
    SPECIES <- list()
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        sppOutSample <- names(samp[i, samp[i, ] < 1])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppOutSample]
            SPECIES[[i]] <- data.frame(sample.dis)
        }
    }
    SPECIES
}


#i = 1
#j = 1
#k = 1
n <- 1000
Site_Rand_Output <- data.frame()
for(i in 1:length(SITENAMES[,1])){
    ii=i
    SITE_OUTPUT <- data.frame()
    tic("Site")
    #SITE_OUTPUT <- data.frame()
    CURRENT_SITE <- SITENAMES[i,]
    print(CURRENT_SITE)
    
    ### Prune Tree to Site
    ABUND_SITE_ALL <- rbind(colnames(ABUND_ALL_SITES),ABUND_ALL_SITES[CURRENT_SITE,])
    SPECIES_SITE <- ABUND_SITE_ALL[1,][which(ABUND_SITE_ALL[2,]>=1)]
    TREE_SITE <- drop.tip(TREE_Obs,which(TREE_Obs$tip.label%in%SPECIES_SITE==FALSE))
    
    ### Prune Abundance to Site
    ABUND_SITE <- ABUND_SITE_ALL[,which(colnames(ABUND_SITE_ALL)%in%TREE_SITE$tip.label)]
    
    ### Stats for Summary Tree
    MNTD_OBSERVED <- Calc.MNTD(ABUND_SITE,IBinomial$IBinomial,TREE = TREE_SITE,n_runs = 999)
    # MNTD_OBSERVED_Between_ALL <- MNTD_OBSERVED[1]
    # MNTD_OBSERVED_Between_NATIVE <- MNTD_OBSERVED[2]
    # MNTD_OBSERVED_Between_NONNATIVE <- MNTD_OBSERVED[3]
    MNTD_NN_N <- MNTD_OBSERVED[4]                                                   #***export*****
    colnames(MNTD_NN_N) <- "MNTD_NN_NPAIR"
    # MNTD_MOD_OBSERVED_Between_N_NN <- MNTD_OBSERVED[5]
    
    ### n Invasive
    I_COUNT <- length(which(TREE_SITE$tip.label%in%IBinomial$IBinomial==TRUE))
    tree_names_n_only <- ABUND_SITE[1,which(colnames(ABUND_SITE)%in%IBinomial$IBinomial==FALSE)]
    
    ### Find N_pair
    temp_out <- Status.0.1.Builder(ABUND_SITE,IBinomial$IBinomial)
    COMM <- temp_out[[1]]
    COMM_ALL <- rbind(COMM,matrix(data=1,nrow = 1,ncol=(length(COMM)/2)))
    #samp <- COMM_ALL
    #dis <- cophenetic(TREE_SITE)
    SPECIES_MATRIX <- MNTD_MOD_return_species_names(COMM_ALL, cophenetic(TREE_SITE), abundance.weighted=FALSE)
    NN_N <- SPECIES_MATRIX[[1]]
    #N_NN <- SPECIES_MATRIX[[2]]
    N_PAIR_LIST <- matrix(NA,I_COUNT,3)
    for (j in 1:I_COUNT){
        ROW <- NN_N[j,]
        N_PAIR_LIST[j,1] <- rownames(ROW)
        N_PAIR_LIST[j,2] <- names(ROW)[which.min(apply(ROW,MARGIN=2,min))]
        N_PAIR_LIST[j,3] <- min(NN_N[j,])
    }
    N_PAIR_LIST <- data.frame(N_PAIR_LIST)
    colnames(N_PAIR_LIST) <- c("NN","N_Pair","Distance")                                         #****export*******
    
    #### Remove NNs from ABUND_SITE and TREE_Site
    ABUND_SITE2 <- ABUND_SITE[,which(colnames(ABUND_SITE)%in%rownames(NN_N)==FALSE)]
    TREE_SITE2 <- drop.tip(TREE_SITE,rownames(NN_N))
    
    ### N_pair to nearest Native
    temp_out_pair <- Status.0.1.Builder(ABUND_SITE2,N_PAIR_LIST$N_Pair)
    COMM_pair <- temp_out_pair[[1]]
    COMM_ALL_pair <- rbind(COMM_pair,matrix(data=1,nrow = 1,ncol=(length(COMM_pair)/2)))
    MNTD_MOD_OUT <- SES_MNTD_MOD(COMM_ALL_pair, cophenetic(TREE_SITE2), abundance.weighted=FALSE)
    #Parse SES.MNTD.MOD
    MNTD_Npair_N <- data.frame(MNTD_MOD_OUT[2,2])    #Compare this to MNTD_NN_N <- MNTD_OBSERVED[4] for site comparisons ***export*******
    colnames(MNTD_Npair_N) <- "MNTD_NPAIR_N"
    
    # Get Distances from Npair to nearest other native
    SPECIES_MATRIX_N_PAIR <- MNTD_MOD_return_species_names(COMM_ALL_pair, cophenetic(TREE_SITE2), abundance.weighted=FALSE)
    NPAIR_N <- SPECIES_MATRIX_N_PAIR[[1]]
    NPAIR_to_N_LIST <- matrix(NA,length(NPAIR_N[,1]),3)
    for (j in 1:length(NPAIR_N[,1])){
        ROW <- NPAIR_N[j,]
        NPAIR_to_N_LIST[j,1] <- rownames(ROW)
        NPAIR_to_N_LIST[j,2] <- names(ROW)[which.min(apply(ROW,MARGIN=2,min))]
        NPAIR_to_N_LIST[j,3] <- min(NPAIR_N[j,])
    }
    NPAIR_to_N_LIST <- data.frame(NPAIR_to_N_LIST)
    colnames(NPAIR_to_N_LIST) <- c("N_Pair","NPair_Nearest_N","Distance")                                     #*****export******
    
    L <- list(MNTD_NN_N,MNTD_Npair_N,N_PAIR_LIST)
    OUT <- do.call(rbind.fill, L)
    L <- list(OUT,NPAIR_to_N_LIST)
    OUT <- do.call(rbind.fill, L)
    
    ### Count
    OBS_ROW <- data.frame()
    for(i in 1:length(N_PAIR_LIST$Distance)){
        ToMatch <- as.character(N_PAIR_LIST[i,2])
        Dist_NN_to_NPair <- as.numeric(as.character(N_PAIR_LIST[i,3]))
        Dist_NPair_to_N <- as.numeric(as.character(NPAIR_to_N_LIST[which(NPAIR_to_N_LIST$N_Pair==ToMatch),3]))
            
        DIFF <- Dist_NN_to_NPair - Dist_NPair_to_N
        OBS_ROW <- rbind(OBS_ROW,DIFF)
            
    }
    NN_FartherFrom_Npair_COUNT <- sum(OBS_ROW >= 0)
    NN_CloserTo_Npair_COUNT <- sum(OBS_ROW < 0)
    if(NN_FartherFrom_Npair_COUNT==NN_CloserTo_Npair_COUNT){RESULT <- "Same"}
    if(NN_FartherFrom_Npair_COUNT > NN_CloserTo_Npair_COUNT){RESULT <- "NN_Farther_From_NPair"}
    if(NN_FartherFrom_Npair_COUNT < NN_CloserTo_Npair_COUNT){RESULT <- "NN_Closer_To_NPair"}
        
    OBS_ROW_Report <- cbind(RESULT,NN_FartherFrom_Npair_COUNT,NN_CloserTo_Npair_COUNT)
    
    ### Write OUT
    filename <- paste("./",DIRECTORY,CURRENT_SITE,"_MNTD_Pair_Analysis_7-7-18.csv",sep="")
    write.csv(OUT,filename,row.names=FALSE)
    
    ######################
    ### Randomizations ###
    ######################
    i=ii
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
    TREE_SITE_wRSP <- drop.tip(TREE_Obs,which(TREE_Obs$tip.label%in%colnames(RSP_AND_SITE_formatted)==FALSE))
    TREE_SITE <- drop.tip(TREE_Obs,which(TREE_Obs$tip.label%in%SPECIES_SITE==FALSE))
    
    # Remove RSP Species that are not in the tree
    RSP_AND_SITE_formatted <- RSP_AND_SITE_formatted[,which(colnames(RSP_AND_SITE_formatted)%in%TREE_SITE_wRSP$tip.label==TRUE)]
    
    ### Prune Abundance to Site
    ABUND_SITE <- ABUND_SITE_ALL[,which(colnames(ABUND_SITE_ALL)%in%TREE_SITE$tip.label)]
    
    ### Stats for Summary Tree
    MNTD_OBSERVED_wRSP <- Calc.MNTD(RSP_AND_SITE_formatted,IBinomial$IBinomial,TREE = TREE_SITE_wRSP,n_runs = 999)
    
    ### Stats for Summary Tree
    MNTD_OBSERVED <- Calc.MNTD(ABUND_SITE,IBinomial$IBinomial,TREE = TREE_SITE,n_runs = 999)
    # MNTD_OBSERVED_Between_ALL <- MNTD_OBSERVED[1]
    # MNTD_OBSERVED_Between_NATIVE <- MNTD_OBSERVED[2]
    # MNTD_OBSERVED_Between_NONNATIVE <- MNTD_OBSERVED[3]
    # MNTD_MOD_OBSERVED_Between_NN-N <- MNTD_OBSERVED[4]
    # MNTD_MOD_OBSERVED_Between_N_NN <- MNTD_OBSERVED[5]
    
    ### n Invasive
    I_COUNT <- length(which(TREE_SITE$tip.label%in%IBinomial$IBinomial==TRUE))
    tree_names_n_only <- ABUND_SITE[1,which(colnames(ABUND_SITE)%in%IBinomial$IBinomial==FALSE)]
    
    MNTD_RAND_Output <- data.frame()
    Distance_Report <- data.frame()
    for (i in 1:n){
        print(paste(i," // ",n,sep=""))
        SIMTREE <- Build.Simulation(TREE_SITE_wRSP,I_COUNT,RSP_SITE_NN_formatted,tree_names_n_only)
        ABUND_SIM <- RSP_AND_SITE_formatted[,colnames(RSP_AND_SITE_formatted)%in%SIMTREE$tip.label==TRUE]
      
        ### Find N_pair
        temp_out <- Status.0.1.Builder(ABUND_SIM,IBinomial$IBinomial)
        COMM <- temp_out[[1]]
        COMM_ALL <- rbind(COMM,matrix(data=1,nrow = 1,ncol=(length(COMM)/2)))
        #samp <- COMM_ALL
        #dis <- cophenetic(SIMTREE)
        SPECIES_MATRIX <- MNTD_MOD_return_species_names(COMM_ALL, cophenetic(SIMTREE), abundance.weighted=FALSE)
        NN_N <- SPECIES_MATRIX[[1]]
        #N_NN <- SPECIES_MATRIX[[2]]
        N_PAIR_LIST <- matrix(NA,I_COUNT,3)
        for (j in 1:I_COUNT){
          ROW <- NN_N[j,]
          N_PAIR_LIST[j,1] <- rownames(ROW)
          N_PAIR_LIST[j,2] <- names(ROW)[which.min(apply(ROW,MARGIN=2,min))]
          N_PAIR_LIST[j,3] <- min(NN_N[j,])
        }
        N_PAIR_LIST <- data.frame(N_PAIR_LIST)
        colnames(N_PAIR_LIST) <- c("NN","N_Pair","Distance")                                         #****export*******
        
        #### Remove NNs from ABUND_SIM and TREE_Site
        ABUND_SITE2 <- ABUND_SIM[,which(colnames(ABUND_SIM)%in%rownames(NN_N)==FALSE)]
        TREE_SITE2 <- drop.tip(SIMTREE,rownames(NN_N))
        
        ### N_pair to nearest Native
        temp_out_pair <- Status.0.1.Builder(ABUND_SITE2,N_PAIR_LIST$N_Pair)
        COMM_pair <- temp_out_pair[[1]]
        COMM_ALL_pair <- rbind(COMM_pair,matrix(data=1,nrow = 1,ncol=(length(COMM_pair)/2)))
        MNTD_MOD_OUT <- SES_MNTD_MOD(COMM_ALL_pair, cophenetic(TREE_SITE2), abundance.weighted=FALSE)
        #Parse SES.MNTD.MOD
        MNTD_Npair_N <- data.frame(MNTD_MOD_OUT[2,2])    #Compare this to MNTD_NN_N <- MNTD_OBSERVED[4] for site comparisons ***export*******
        colnames(MNTD_Npair_N) <- "MNTD_NPAIR_N"
        
        # Get Distances from Npair to nearest other native
        SPECIES_MATRIX_N_PAIR <- MNTD_MOD_return_species_names(COMM_ALL_pair, cophenetic(TREE_SITE2), abundance.weighted=FALSE)
        NPAIR_N <- SPECIES_MATRIX_N_PAIR[[1]]
        NPAIR_to_N_LIST <- matrix(NA,length(NPAIR_N[,1]),3)
        for (j in 1:length(NPAIR_N[,1])){
          ROW <- NPAIR_N[j,]
          NPAIR_to_N_LIST[j,1] <- rownames(ROW)
          NPAIR_to_N_LIST[j,2] <- names(ROW)[which.min(apply(ROW,MARGIN=2,min))]
          NPAIR_to_N_LIST[j,3] <- min(NPAIR_N[j,])
        }
        NPAIR_to_N_LIST <- data.frame(NPAIR_to_N_LIST)
        colnames(NPAIR_to_N_LIST) <- c("N_Pair","NPair_Nearest_N","Distance")                                     #*****export******
        
        # L <- list(MNTD_NN_N,MNTD_Npair_N,N_PAIR_LIST)
        # OUT <- do.call(rbind.fill, L)
        # L <- list(OUT,NPAIR_to_N_LIST)
        # OUT <- do.call(rbind.fill, L)
        
        RAND_ROW <- data.frame()
        for(i in 1:length(N_PAIR_LIST$Distance)){
            ToMatch <- as.character(N_PAIR_LIST[i,2])
            Dist_NN_to_NPair <- as.numeric(as.character(N_PAIR_LIST[i,3]))
            Dist_NPair_to_N <- as.numeric(as.character(NPAIR_to_N_LIST[which(NPAIR_to_N_LIST$N_Pair==ToMatch),3]))
            
            DIFF <- Dist_NN_to_NPair - Dist_NPair_to_N
            RAND_ROW <- rbind(RAND_ROW,DIFF)
        }
        NN_FartherFrom_Npair_COUNT <- sum(RAND_ROW >= 0)
        NN_CloserTo_Npair_COUNT <- sum(RAND_ROW < 0)
        if(NN_FartherFrom_Npair_COUNT==NN_CloserTo_Npair_COUNT){RESULT <- "Same"}
        if(NN_FartherFrom_Npair_COUNT > NN_CloserTo_Npair_COUNT){RESULT <- "NN_Farther_From_NPair"}
        if(NN_FartherFrom_Npair_COUNT < NN_CloserTo_Npair_COUNT){RESULT <- "NN_Closer_To_NPair"}
        
        RAND_ROW_Report <- cbind(RESULT,NN_FartherFrom_Npair_COUNT,NN_CloserTo_Npair_COUNT)
        Distance_Report <- rbind(Distance_Report,RAND_ROW_Report)
      
        SIMTREE_MNTD <- data.frame(Calc.MNTD(ABUND_SIM,IBinomial$IBinomial,TREE = SIMTREE,n_runs = 1))
        colnames(SIMTREE_MNTD) <- c("MNTD_ALL_species_at_site",
                                    "MNTD_between_NATIVE_species_at_site",
                                    "MNTD_between_NONNATIVE_species_at_site",
                                    "MNTD_MOD_NN_N",
                                    "MNTD_MOD_N_NN")
        MNTD_RAND_Output <- rbind.data.frame(MNTD_RAND_Output,SIMTREE_MNTD)
    }
    
    colnames(Distance_Report) <- c("Result","NN_Far_From_NPair","NN_Close_To_NPair")
    
    SAME_TOTAL <- sum(Distance_Report$Result=="Same")/1000
    NN_FAR_TOTAL <- sum(Distance_Report$Result=="NN_Farther_From_NPair")/1000
    NN_CLOSE_TOTAL <- sum(Distance_Report$Result=="NN_Closer_To_NPair")/1000
    
    if(as.numeric(as.character(OBS_ROW_Report[1,2])) <= I_COUNT/2){
        P_NN_Far <- sum(as.numeric(as.character(Distance_Report[,2])) > as.numeric(as.character(OBS_ROW_Report[1,2])))/1000
        P_NN_Close <- sum(as.numeric(as.character(Distance_Report[,3])) >= as.numeric(as.character(OBS_ROW_Report[1,3])))/1000
    }else{
        P_NN_Far <- sum(as.numeric(as.character(Distance_Report[,2])) <= as.numeric(as.character(OBS_ROW_Report[1,2])))/1000
        P_NN_Close <- sum(as.numeric(as.character(Distance_Report[,3])) < as.numeric(as.character(OBS_ROW_Report[1,3])))/1000
    }
    H <- hist(as.numeric(as.character(Distance_Report[,2])))
    quantile(as.numeric(as.character(Distance_Report[,2])))
    H_Pvalue <- Calc.P.Value(as.numeric(as.character(OBS_ROW_Report[1,2])),as.numeric(as.character(Distance_Report[,2])),1000)
    if(H_Pvalue>.5){H_Pvalue <- 1 - H_Pvalue}
    Site_Rand_Output <- rbind(Site_Rand_Output,cbind(as.character(CURRENT_SITE),OBS_ROW_Report,SAME_TOTAL,NN_FAR_TOTAL,NN_CLOSE_TOTAL,H_Pvalue))#P_NN_Far,P_NN_Close))

}#END Site Loop
colnames(Site_Rand_Output) <- c("Site","MCC_Result","NN_Farther_From_NPair","NN_Closer_To_NPair",
                                "Equal","NN_Far_From_NPair","NN_Close_To_NPair","Hist_P-Value")#,"P_NN_Far","P_NN_Close")
filename2 <- paste("./",DIRECTORY,"MNTD_Pair_Analysis_Randomizations_P-Values_7-7-18.csv",sep="")
write.csv(Site_Rand_Output,filename2,row.names=FALSE)
