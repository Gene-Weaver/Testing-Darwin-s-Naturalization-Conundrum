# Code for: https://dx.doi.org/10.1111/ddi.12861
#
#     Ng, J., W. Weaver,* and R.G. Laport. In Press. 
#         Testing Darwin's naturalization conundrum using phylogenetic relationships: generalizable 
#         patterns across disparate communities. Diversity and Distributions. 

###### Set the directory location for output files. If it is in the working directory use ""
DIRECTORY <- ""

library(ape)
library(picante)
require(phytools)
library(evobiR)
library(phyloch)

source("./Toolbox/Status.0.1.Builder.R")
source("./Toolbox/Observed.PAE.R")
source("./Toolbox/Observed.PAE.Fast.R")
source("./Toolbox/Build.Simulation.R")
source("./Toolbox/Simulated.PAE.R")
source("./Toolbox/Observed.MNTD.R")
source("./Toolbox/Simulated.MNTD.R")

################################################################
###### Import Abundance Data with Species as Columns ###########
################################################################
# Data
ABUND_RAW <- read.csv("./Abundance_Data_Species_as_Col_SEVEN.csv")
ABUND_ALL_SITES <- ABUND_RAW[,-1]

# Get Site Names
SITENAMES <- data.frame(as.character(ABUND_RAW[,1]))
colnames(SITENAMES) <- "Site"

#####################################################
###### *Import Tree and Invasive Status* ############
#####################################################
# June 2018 tree
TREE_RAW_0 <- read.nexus("./NEONtree_wgymnoclade_1000_NEW.trees")

# Import USDA Invasive csv
INVASIVE_0 <- read.csv("./USDA_Invasive.csv")
###################################################################
######## Invasive Binomial Builder ################################
###################################################################
# Construct Binomial List for non-native species, INVASIVE is the code not the binomial name
temp <- Invasive.Binomial.Builder(INVASIVE_0)
IBinomial <- data.frame(temp$IBinomial)
INVASIVE <- data.frame(temp$Invasive)
colnames(IBinomial) <- c("IBinomial")
colnames(INVASIVE) <- c("Invasive")

###Build list of NEON invasive species for random sampling in MNTD
ALL_NEON_SPECIES <- data.frame(colnames(ABUND_ALL_SITES))
colnames(ALL_NEON_SPECIES) <- c("ALL_NEON_SPECIES")
INVASIVE_NEON_SPECIES <- data.frame(ALL_NEON_SPECIES$ALL_NEON_SPECIES[ALL_NEON_SPECIES$ALL_NEON_SPECIES %in% IBinomial$IBinomial])
colnames(INVASIVE_NEON_SPECIES) <- c("INVASIVE_NEON_SPECIES")

#TOTAL_SPECIES_LIST <- data.frame(TREE_RAW_0[[1]]$tip.label[which(TREE_RAW_0[[1]]$tip.label%in%ALL_NEON_SPECIES$ALL_NEON_SPECIES)])

# ***Loop through each site***
for(i in 1:length(SITENAMES[,1])){
    # Define Site Name
    CURRENT_SITE <- SITENAMES[i,]
    print(CURRENT_SITE)
    # Track Progress
    j <- 1
    tree_index <- 1
    # abund = site-specific abundance data 
    abund <- rbind(colnames(ABUND_ALL_SITES),ABUND_ALL_SITES[CURRENT_SITE,])
    
    # Initialize Storage Dataframes for Output Data
    # PAE
    EXPORT_SESPD_OBSERVED <- data.frame()
    EXPORT_SESPD_GROUPS <- data.frame()
    EXPORT_PAE_ALL <- data.frame()
    TREE_P_VALUE <- data.frame()
    # MNTD
    EXPORT_MNTD_OBSERVED <- data.frame()
    TREE_P_VALUE_MNTD <- data.frame()
    # SESPD
    TREE_P_VALUE_PD <- data.frame()
    
    #######################################
    ##### Define Site-Specific Species #### ***This assumes that each of the 1,000 trees has identical species***
    ####################################### all _DEFINE variable should only be used for species names
    tree_DEFINE <- TREE_RAW_0[[1]]
    
    # Names of all species in tree
    Tree_Names_DEFINE <- tree_DEFINE$tip.label
    # Site-Specific Species based on having an abundance >0
    Site_Names_DEFINE <- abund[1,][which(abund[2,]>=1)]
    Site_Size_DEFINE <- length(Site_Names_DEFINE)
    
    # Groups of Species
    REMOVE_THESE <- Tree_Names_DEFINE[which((Tree_Names_DEFINE%in%Site_Names_DEFINE)==FALSE)]
    TREE_NAMES_SITE_SPECIFIC <- Tree_Names_DEFINE[which(Tree_Names_DEFINE%in%Site_Names_DEFINE)]
    SPECIES_NOT_IN_TREE <- Site_Names_DEFINE[which((Site_Names_DEFINE%in%TREE_NAMES_SITE_SPECIFIC)==FALSE)]
    # Define Non-native species 
    REMOVE_I <- as.character(colnames(Site_Names_DEFINE[Site_Names_DEFINE%in%IBinomial$IBinomial]))
    
    # ***Loop through each Tree*** there are 1,000 trees in final SEVEN analysis
    for(k in 1: length(TREE_RAW_0)){
        tree_cycle <- TREE_RAW_0[[k]]
        # Track Progress
        print(paste(CURRENT_SITE,tree_index),sep="  ")
        
        #############################################
        ###### Drop Tips - Site Specific ############
        #############################################
        tree <- tree_cycle
        
        # Site-specific tree is TREE, while tree is the whole unmodified tree with 319ish species
        TREE <- drop.tip(tree,REMOVE_THESE)
        # Maintain copy of unmodified observed tree
        TREE_0 <- TREE
        
        #############################
        ###### Drop Invasives #######
        #############################    
        TREE_NAMES <- TREE$tip.label
        TREE_N_ONLY <- drop.tip(TREE,REMOVE_I)
        # Maintain unmodified native only tree
        TREE_N_ONLY_O <- TREE_N_ONLY
        tree_names_n_only <- TREE_N_ONLY_O$tip.label
        
        # abund site-specific based on species present in the tree
        ABUND <- abund[,which(colnames(abund)%in%TREE_NAMES)]
        ABUND_N_ONLY <- abund[,which(colnames(abund)%in%tree_names_n_only)]
        
        # make 0/1 invasive/native for sespd
        #Native Status Builder
        temp_out <- Status.0.1.Builder(ABUND,IBinomial$IBinomial)
        COMM <- temp_out[[1]]
        COMM_N_ONLY <- temp_out[[2]]
        
        ######################################
        ###### Calculate Observed PAE ########
        ######################################
        # temp_out[[1]] = EXPORT_SESPD <- cbind(SESPD,SESPD_N_ONLY)  *1 row
        # temp_out[[2]] = EXPORT_SESPD_GROUPS <- SESPD_Groups        *2 rows
        # temp_out[[3]] = EXPORT_PAE_ALL <- cbind(PDvalue_0,PAEdividend,PAEdivisor,PAEquotient,OVERALL_PAE_RATIO)
        # temp_out[[4]] = EXPORT_PAE_N_ONLY <- cbind(PDvalueN_0,PAEdividendN,PAEdivisorN,PAEquotientN)
        # temp_out[[5]] = OVERALL_PAE_RATIO
        temp_out <- Observed.PAE.Fast(TREE,TREE_N_ONLY,COMM,ABUND,ABUND_N_ONLY)
        EXPORT_SESPD_OBSERVED <-  rbind(EXPORT_SESPD_OBSERVED,temp_out[[1]])
        EXPORT_SESPD_GROUPS <- rbind(EXPORT_SESPD_GROUPS,cbind(rbind(paste(CURRENT_SITE,"_NN",sep=""),paste(CURRENT_SITE,"_N",sep="")),temp_out[[2]]))
        EXPORT_PAE_ALL <- rbind(EXPORT_PAE_ALL,cbind(temp_out[[3]],temp_out[[4]])) #overallPAEratio seperates the two
        
        
        ######################################
        ###### Calculate Observed MNTD #######
        ######################################
        # temp_out_MNTD[[1]] --> 2 row export values
        # temp_out_MNTD[[2]] --> SES_MNTD_MOD_ratio
        
        temp_out_MNTD <- Observed.MNTD(TREE,COMM)
        EXPORT_MNTD_OBSERVED <- rbind(EXPORT_MNTD_OBSERVED,cbind(rbind(paste(as.character(CURRENT_SITE),tree_index,sep="_"),paste(as.character(CURRENT_SITE),tree_index,sep="_")),temp_out_MNTD[[1]]))
        OBSERVED_MNTD_RATIO <- temp_out_MNTD[[2]]
        
        ######################################
        #### Simulation for PAE and MNTD #####
        ######################################
        print(paste(as.character(CURRENT_SITE),tree_index,"Simulation",sep=" "))
        numReps <- 2
        rep<-1
        
        ### Define p-value matrices
        #PAE
        SIMULATION_PAE <- matrix(nrow=numReps, ncol=1)
        colnames(SIMULATION_PAE)<-"PAE_Ratio"
        #MNTD
        #SIMULATION_MNTD <- matrix(nrow=numReps, ncol=1)
        #colnames(SIMULATION_MNTD)<-"MNTD_Ratio"
        #SESPD
        SIMULATION_PD <- matrix(nrow=numReps, ncol=1)
        colnames(SIMULATION_PD)<-"PD_Ratio"
        
        SITE_INVASIVE_COUNT <- length(REMOVE_I)
        # *Start Simulation*
        while (rep <= numReps){
            print(rep)
            # temp_sim_out[[1]] --> ABUND_S
            # temp_sim_out[[2]] --> ABUND_N_ONLY_S
            # temp_sim_out[[3]] --> sim_tree_RAND
            # temp_sim_out[[4]] --> sim_tree_RAND_N_ONLY
            temp_sim_out <- Build.Simulation(tree,SITE_INVASIVE_COUNT)
            
            # Run PAE Simulation
            #print(paste("PAE Start",rep,sep=" "))
            SIM_PAE <- Simulated.PAE(temp_sim_out[[3]],temp_sim_out[[4]],temp_sim_out[[1]],temp_sim_out[[2]],IBinomial$IBinomial)
            #print("PAE End")
            # SIM_PAE[[1]] --> OVERALL_PAE_RATIO
            # SIM_PAE[[2]] --> EXPORT_PAE = cbind(SESPD,SESPD_N_ONLY,PDvalue_0,PAEdividend,PAEdivisor,PAEquotient,OVERALL_PAE_RATIO,PDvalueN_0,PAEdividendN,PAEdivisorN,PAEquotientN)
            # SIM_PAE[[3]] --> EXPORT_SESPD_GROUPS
            SIMULATION_PAE[rep,1] <- SIM_PAE[[1]]
            
            # Run MNTD Simulation
            #print(paste("MNTD Start",rep,sep=" "))
            SIM_MNTD <- Simulated.MNTD(temp_sim_out[[3]],temp_sim_out[[1]],IBinomial$IBinomial)
            #print("MNTD End")
            SIMULATION_MNTD[rep,1] <- SIM_MNTD
            
            # PD Simulation for p-value
                # SIM_PAE[[2]][2] is the pd.obs for ses.pd with all species
                # SIM_PAE[[2]][10] is the pd.obs for ses.pd with only native species
            SIMULATION_PD[rep,1] <- as.numeric(SIM_PAE[[2]][2]/SIM_PAE[[2]][10])
            
            rep <- rep+1
        }
        # PAE: Get p-value and 95% CI for null simulation based on PAE Ratio
        p_v <- sum(SIMULATION_PAE > temp_out[[5]])/numReps
        if(p_v >.5){p_v <- 1-p_v}
        d <- density(SIMULATION_PAE[,1])
        q2_5 <- quantile(SIMULATION_PAE[,1], 0.025)
        q97_5 <- quantile(SIMULATION_PAE[,1], 0.975)
        x1 <- min(d$x[which(d$x > q2_5)])
        x2 <- max(d$x[which(d$x < q97_5)])
        
        # MNTD: Get p-value and 95% CI for null simulation based on PAE Ratio
        p_vM <- sum(SIMULATION_MNTD > temp_out_MNTD[[2]])/numReps
        if(p_vM >.5){p_vM <- 1-p_vM}
        dM <- density(SIMULATION_MNTD[,1])
        q2_5M <- quantile(SIMULATION_MNTD[,1], 0.025)
        q97_5M <- quantile(SIMULATION_MNTD[,1], 0.975)
        x1M <- min(dM$x[which(dM$x > q2_5M)])
        x2M <- max(dM$x[which(dM$x < q97_5M)])
        
        # PD: Get p-value and 95% CI for null simulation based on PAE Ratio
        SESPD_RATIO_OBSERVED <- temp_out[[1]][2]/temp_out[[1]][10]
        p_vPD <- sum(SIMULATION_PD > SESPD_RATIO_OBSERVED)/numReps
        if(p_vPD >.5){p_vPD <- 1-p_vPD}
        dP <- density(SIMULATION_PD[,1])
        q2_5P <- quantile(SIMULATION_PD[,1], 0.025)
        q97_5P <- quantile(SIMULATION_PD[,1], 0.975)
        x1P <- min(dP$x[which(dP$x > q2_5P)])
        x2P <- max(dP$x[which(dP$x < q97_5P)])
        
        # Determine if simulation was significant
        if(temp_out[[5]] > x1 & temp_out[[5]] < x2){PAE_SIG <- NA}else{PAE_SIG <- "SIGNIFICANT"}
        if(temp_out_MNTD[[2]] > x1M & temp_out_MNTD[[2]] < x2M){MNTD_SIG <- NA}else{MNTD_SIG <- "SIGNIFICANT"}
        if(SESPD_RATIO_OBSERVED > x1P & SESPD_RATIO_OBSERVED < x2P){PD_SIG <- NA}else{PD_SIG <- "SIGNIFICANT"}
        
        # Current site,tree index, observed PAE,p-value, lower CI, upper CI
        TREE_P_VALUE <- rbind(TREE_P_VALUE,cbind(as.character(CURRENT_SITE),tree_index,as.character(PAE_SIG),temp_out[[5]],p_v,x1,x2))
        TREE_P_VALUE_MNTD <- rbind(TREE_P_VALUE_MNTD,cbind(as.character(CURRENT_SITE),tree_index,as.character(MNTD_SIG),temp_out_MNTD[[2]],p_vM,x1M,x2M))
        TREE_P_VALUE_PD <- rbind(TREE_P_VALUE_PD,cbind(as.character(CURRENT_SITE),tree_index,as.character(PD_SIG),SESPD_RATIO_OBSERVED,p_vPD,x1P,x2P))
        
        tree_index <- tree_index + 1
    }#END TREE LOOP
    
    # Colnames p-values
    colnames(TREE_P_VALUE) <- c("Site","Tree_Iteration","Significance","Observed_PAE","p_value","LowerCI","UpperCI")
    colnames(TREE_P_VALUE_MNTD) <- c("Site","Tree_Iteration","Significance","Observed_PAE","p_value","LowerCI","UpperCI")
    colnames(TREE_P_VALUE_PD) <- c("Site","Tree_Iteration","Significance","Observed_SESPD","p_value","LowerCI","UpperCI")
    
    # Filenames PAE
    filename_TREE_P_VALUE <- paste("./",DIRECTORY,CURRENT_SITE,"_PAE_Simulation_PValues.csv",sep="")
    filename_EXPORT_SESPD_OBSERVED <- paste("./",DIRECTORY,CURRENT_SITE,"_SESPD_Observed.csv",sep="")
    filename_EXPORT_SESPD_GROUPS <- paste("./",DIRECTORY,CURRENT_SITE,"_SESPD_Groups.csv",sep="")
    filename_EXPORT_PAE_ALL <- paste("./",DIRECTORY,CURRENT_SITE,"_PAE_Observed.csv",sep="")
    
    # Filename PD p-values
    filename_TREE_P_VALUE_PD <- paste("./",DIRECTORY,CURRENT_SITE,"_SESPD_Simulation_PValues.csv",sep="")
    write.csv(TREE_P_VALUE_PD,filename_TREE_P_VALUE_PD,row.names=FALSE)
    
    # Write PAE: TREE_P_VALUE, EXPORT_SESPD_OBSERVED, EXPORT_SESPD_GROUPS, EXPORT_PAE_ALL  
    write.csv(TREE_P_VALUE,filename_TREE_P_VALUE,row.names=FALSE)
    write.csv(EXPORT_SESPD_OBSERVED,filename_EXPORT_SESPD_OBSERVED,row.names=FALSE)
    write.csv(EXPORT_SESPD_GROUPS,filename_EXPORT_SESPD_GROUPS,row.names=FALSE)
    write.csv(EXPORT_PAE_ALL,filename_EXPORT_PAE_ALL,row.names=FALSE)
    
    # Colnames MNTD
    colnames(EXPORT_MNTD_OBSERVED) <- c("Tree_Iteration","MNTD reg. (A)","MNTD MOD (B)","MNTD ratio (C)","MNTD MOD ratio (D)","SES MNTD MOD ratio (E)","SES MNTD reg (I)","SES MNTD MOD (J)","SES MOD sd (G)","SES reg sd (H)","SES_sd_ratio [SES_MOD_sd/SES_REG_sd] (F)","START SES MOD n_taxa","mntd.obs","mntd.rand.mean","mntd.rand.sd","mntd.obs.rank","mntd.obs.z","mntd.obs.p","runs",
                                        "START SES reg n_taxa","mntd.obs","mntd.rand.mean","mntd.rand.sd","mntd.obs.rank","mntd.obs.z","mntd.obs.p","runs")
    
    # Filenames MNTD 
    filename_TREE_P_VALUE_MNTD <- paste("./",DIRECTORY,CURRENT_SITE,"_MNTD_Simulation_PValues.csv",sep="")
    filename_EXPORT_MNTD_OBSERVED <- paste("./",DIRECTORY,CURRENT_SITE,"_MNTD_Observed.csv",sep="")
    # Write MNTD
    write.csv(TREE_P_VALUE_MNTD,filename_TREE_P_VALUE_MNTD,row.names=FALSE)
    write.csv(EXPORT_MNTD_OBSERVED,filename_EXPORT_MNTD_OBSERVED,row.names=FALSE)
    
}#END SITE LOOP
