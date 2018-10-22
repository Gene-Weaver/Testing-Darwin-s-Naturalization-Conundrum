# Code for: https://dx.doi.org/10.1111/ddi.12861
#
#     Ng, J., W. Weaver,* and R.G. Laport. In Press. 
#         Testing Darwin's naturalization conundrum using phylogenetic relationships: generalizable 
#         patterns across disparate communities. Diversity and Distributions. 

library(ape)
library(picante)
require(phytools)
library(evobiR)
library(phyloch)
source("./Toolbox/Invasive_Binomial_Builder.R")
source("./Toolbox/Prune_Tree.R")

###### Set the directory location for output files. If it is in the working directory use ""
DIRECTORY <- ""

################################################################
###### Phylogenetic Abundance Evenness PAE - ALL SITES #########
################################################################
#           PAE Equation From: http://rfunctions.blogspot.com/2013/07/pae-function-phylogenetic-abundance.html
#           PAE falls between 0 and 1 when individuals are clustered into relatively short terminal branches, 
#           whereas PAE > 1 when individuals are clustered into long terminal branches. 
#           When individuals are evenly distributed across terminal branch lengths, PAE = 1.

################################################################
###### Import Abundance Data with Species as Columns ###########
################################################################
#ABUND_RAW <- read.csv("./FOR_TESTING2_Abundance_Data_Species_as_Col.csv")
ABUND_RAW <- read.csv("./Abundance_Data_Species_as_Col_SEVEN.csv")
ABUND_ALL_SITES <- ABUND_RAW[,-1]

SITENAMES <- data.frame(as.character(ABUND_RAW[,1]))
colnames(SITENAMES) <- "Site"
################################################################
###### Import Tree, Drop Tips Based on Site Species ############
################################################################
TREE_RAW_0 <- read.nexus("./NEONtree_wgymnoclade_1000_NEW.trees")

### List of species that should not appear in randomizations + one species that was duplicated
IGNORE <- data.frame(c("Macfadyena_unguis.cati","Smilax_bona.nox","Abies_alba",
                       "Pinus_thunbergii","Callitris_glaucophylla","Platycladus_orientalis",
                       "Taxus_baccata","Taxus_cuspidata","Lygodium_microphyllum",
                       "Lygodium_japonicum"))

INVASIVE_0 <- read.csv("./USDA_Invasive.csv")
###################################################################
######## Invasive Binomial Builder ################################
###################################################################
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

################################################################
###### Calculate PAE  ##########################################
################################################################

for(i in 1:length(SITENAMES[,1])){    
    CURRENT_SITE <- SITENAMES[i,]
    print(CURRENT_SITE)
    abund <- rbind(colnames(ABUND_ALL_SITES),ABUND_ALL_SITES[i,])
    EXPORT_PAE <- data.frame()
    EXPORT_SESPD <- data.frame()
    EXPORT_SESPD_Groups <- data.frame()

    tree_index <- 1
    EXPORT_ALL_PAE <- data.frame()
    EXPORT_ALL_SESPD <- data.frame()
    for(k in 1:length(TREE_RAW_0)){
        tree_cycle <- TREE_RAW_0[[k]]
      
        print(paste(CURRENT_SITE,tree_index),sep="  ")
        tree_index <- tree_index + 1

        #############################################
        ###### Drop Tips - Site Specific ############
        #############################################
        tree <- tree_cycle
        
        Tree_Names <- tree$tip.label
        # Site-specific Names Based on NEON Data, are refined later by removing species that are not in the tree
        Site_Names <- abund[1,][which(abund[2,]>=1)]
        Site_Size <- length(Site_Names)
        
        #Tree_Size <- length(tree$tip.label)
        REMOVE_THESE <- Tree_Names[which((Tree_Names%in%Site_Names)==FALSE)]
        SPECIES_NOT_IN_TREE <- Site_Names[which((Site_Names%in%REMOVE_THESE)==TRUE)]
        TREE_NAMES_SITE_SPECIFIC <- Tree_Names[which(Tree_Names%in%Site_Names)]
        
        # Site-specific tree is TREE, while tree is the whole unmodified tree with 319ish species
        TREE <- drop.tip(tree,REMOVE_THESE)
        TREE <- drop.tip(TREE,as.character(IGNORE[,1]))
        TREE_0 <- TREE
        
        #############################
        ###### Drop Invasives #######
        #############################    
        TREE_NAMES <- TREE$tip.label
        REMOVE_I <- TREE_NAMES[TREE_NAMES%in%IBinomial$IBinomial]
        TREE_N_ONLY <- drop.tip(TREE,REMOVE_I)
        TREE_N_ONLY_O <- TREE_N_ONLY
        tree_names_n_only <- TREE_N_ONLY_O$tip.label
        
        # abund site specific
        ABUND <- abund[,which(colnames(abund)%in%TREE_NAMES)]
        ABUND_N_ONLY <- abund[,which(colnames(abund)%in%tree_names_n_only)]
        # make 0/1 invasive/native for sespd
        # Native Status Builder
        tree_names_STAT_1 <- colnames(ABUND)
        tree_names_STAT_2 <- colnames(ABUND)
        tree_names_STAT_1[tree_names_STAT_1 %in% IBinomial$IBinomial] <- 1
        tree_names_STAT_1[(tree_names_STAT_1 ==1)==FALSE] <- 0
        tree_names_STAT_2[tree_names_STAT_2 %in% IBinomial$IBinomial] <- 0
        tree_names_STAT_2[(tree_names_STAT_2 ==0)==FALSE] <- 1
        COMM <- rbind(tree_names_STAT_1,tree_names_STAT_2)
        colnames(COMM) <- colnames(ABUND)
        COMM_N_ONLY <- COMM[,(which(colnames(COMM)%in%IBinomial$IBinomial==FALSE))]
        
        #############################
        ###### Calculate PAE ########
        #############################
        terms <- TREE$edge[,2] <= Ntip(TREE)
        terminal.edges <- TREE$edge.length[terms]
        names(terminal.edges)<-TREE$tip.label
        PDvalue<-sum(TREE$edge.length)
        PDvalue_0 <- PDvalue
        
        SESPD_Groups <- ses.pd(COMM,TREE,null.model = "taxa.labels", runs = 999, iterations = 1000)
        SESPD <- ses.pd(ABUND,TREE,null.model = "taxa.labels", runs = 999, iterations = 1000)
        SESPD <- SESPD[1,]
        
        spnames<-names(terminal.edges)
        termabun<-numeric(length(terminal.edges))
        for(i in 1:length(terminal.edges)){
            sp<-spnames[i] # get species names
            termabun[i]<-(terminal.edges[sp])*(as.integer((ABUND[sp][2,]))-1)
        }
        PAEdividend<-(sum(unlist(termabun)))+PDvalue
        #Original
        #PAEdivisor<-PDvalue+(((sum(as.integer(abund[2,][which(abund[2,]>=0)]))/length(abund[which(abund[2,]>=0)]))-1)*sum(terminal.edges))
        PAEdivisor<-PDvalue+(((sum(as.integer(ABUND[2,][which(ABUND[2,]>=1)]))/length(ABUND[which(ABUND[2,]>=1)]))-1)*sum(terminal.edges))
        PAEquotient<-PAEdividend/PAEdivisor
        print(PAEquotient)
        
        ############################
        ## Calculate PAE - N Only ##
        ############################
        termsN <- TREE_N_ONLY$edge[,2] <= Ntip(TREE_N_ONLY)
        terminal.edgesN <- TREE_N_ONLY$edge.length[termsN]
        names(terminal.edgesN)<-TREE_N_ONLY$tip.label
        PDvalueN<-sum(TREE_N_ONLY$edge.length)
        PDvalueN_0 <- PDvalueN
        
        SESPD_N_ONLY <- ses.pd(ABUND_N_ONLY,TREE_N_ONLY,null.model = "taxa.labels", runs = 999, iterations = 1000)
        SESPD_N_ONLY <- SESPD_N_ONLY[1,]
        
        spnamesN<-names(terminal.edgesN)
        termabunN<-numeric(length(terminal.edgesN))
        for(i in 1:length(terminal.edgesN)){
            spN<-spnamesN[i] # get species names
            termabunN[i]<-(terminal.edgesN[spN])*(as.integer((ABUND_N_ONLY[spN][2,]))-1)
        }
        PAEdividendN<-(sum(unlist(termabunN)))+PDvalueN
        PAEdivisorN<-PDvalueN+(((sum(as.integer(ABUND_N_ONLY[2,][which(ABUND_N_ONLY[2,]>=1)]))/length(ABUND_N_ONLY[which(ABUND_N_ONLY[2,]>=1)]))-1)*sum(terminal.edgesN))
        PAEquotientN<-PAEdividendN/PAEdivisorN
        print(PAEquotientN)
    
    
        OVERALL_PAE_RATIO <- PAEquotient/PAEquotientN
        PAEquotient_O <- PAEquotient
        PAEquotientN_O <- PAEquotientN
    
        ################################################################################################################################
        ###### Shuffle tips for null distribution SIMULATION ###########################################################################
        ################################################################################################################################
        numReps <- 1000
        rep<-1
        SIMULATION <- matrix(nrow=numReps, ncol=1)
        #SIMULATION_SESPD <- matrix(nrow=numReps, ncol=1)
        colnames(SIMULATION)<-"PAE_Ratio"
        #colnames(SIMULATION_SESPD)<-"SESPD"
        SITE_INVASIVE_COUNT <- length(REMOVE_I)
        while (rep <= numReps){
            #BUILD SIM TREE
            SIM_TREE <- tree
            # take random sample of invasive species and add the names to sim_species 
            sim_I_SAMPLE_0 <- data.frame(sample(INVASIVE_NEON_SPECIES$INVASIVE_NEON_SPECIES,SITE_INVASIVE_COUNT,replace = FALSE))
            #NAtive names carried through simulations per site
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
            colnames(sim_N) <- "Sim_Species_All"
            #merge sim_I and natives
            sim_species <- rbind.data.frame(sim_I_SAMPLE[1],sim_N[1])
            #### drop all other tips to make it site specific with the random Invasives 
            #setdiff(tree$tip.label, species)
            sim_tree_RAND <- drop.tip(SIM_TREE,setdiff(SIM_TREE$tip.label, sim_species$Sim_Species_All))
            #sim_tree_RAND <- drop.tip(SIM_TREE,which(SIM_TREE$tip.label%in%sim_species$Sim_Species_All==FALSE))
            sim_tree_RAND_N_ONLY <- drop.tip(sim_tree_RAND,which(sim_tree_RAND$tip.label%in%INVASIVE_NEON_SPECIES$INVASIVE_NEON_SPECIES))
            
            #Make abund for each sim
            ABUND_S <- abund[,which(colnames(abund)%in%sim_tree_RAND$tip.label)]
            ABUND_N_ONLY_S <- abund[,which(colnames(abund)%in%sim_tree_RAND_N_ONLY$tip.label)]
            
            #### CHECK lengths ######
            print(paste(CURRENT_SITE,"-",tree_index,"  ",round((rep/numReps)*100,2),"% total number of species ",length(sim_species$Sim_Species_All),sep=''))
            print(paste("total number of Invasive species ",length(sim_I_SAMPLE$Sim_Species_All),sep=''))
            print(paste("rand tree ",length(sim_tree_RAND$tip.label),sep=""))
            print(paste("rand n only ",length(sim_tree_RAND_N_ONLY$tip.label),sep=""))
            #############################
            #### Calculate PAE Sim ######
            #############################
            #TREE <- tipShuffle(TREE)
            #sim_tree_rand
            terms <- sim_tree_RAND$edge[,2] <= Ntip(sim_tree_RAND)
            terminal.edges <- sim_tree_RAND$edge.length[terms]
            names(terminal.edges)<-sim_tree_RAND$tip.label
            
            # SESPD Sim 
            #SESPD_Sim <- ses.pd(ABUND_S,sim_tree_RAND,null.model = "taxa.labels", runs = 999, iterations = 1000)
            #SESPD_Sim <- SESPD_Sim[1,]
            #SESPD_Sim_Value <- SESPD_Sim[1,2]
            
            # Total phylogenetic distance for calculating PAE
            PDvalue<-sum(sim_tree_RAND$edge.length)
            
            spnames<-names(terminal.edges)
            termabun<-numeric(length(terminal.edges))
            for(i in 1:length(terminal.edges)){
                sp<-spnames[i] # get species names
                termabun[i]<-(terminal.edges[sp])*(as.integer((ABUND_S[sp][2,]))-1)
            }
            PAEdividend<-(sum(unlist(termabun)))+PDvalue
            PAEdivisor<-PDvalue+(((sum(as.integer(ABUND_S[2,][which(ABUND_S[2,]>=1)]))/length(ABUND_S[which(ABUND_S[2,]>=1)]))-1)*sum(terminal.edges))
            PAEquotient <- PAEdividend/PAEdivisor
    
            ###########################################
            ####### Calculate PAE Sim - N Only ########
            ###########################################
            termsN <- sim_tree_RAND_N_ONLY$edge[,2] <= Ntip(sim_tree_RAND_N_ONLY)
            terminal.edgesN <- sim_tree_RAND_N_ONLY$edge.length[termsN]
            names(terminal.edgesN)<-sim_tree_RAND_N_ONLY$tip.label
            
            #SESPD_N_ONLY_Sim <- ses.pd(ABUND_N_ONLY_S,sim_tree_RAND_N_ONLY,null.model = "taxa.labels", runs = 999, iterations = 1000)
            #SESPD_N_ONLY_Sim <- SESPD_N_ONLY_Sim[1,]
            #SESPD_N_ONLY_Sim_Value <- SESPD_N_ONLY_Sim[1,2]
            
            PDvalueN<-sum(sim_tree_RAND_N_ONLY$edge.length)
            
            spnamesN<-names(terminal.edgesN)
            termabunN<-numeric(length(terminal.edgesN))
            for(i in 1:length(terminal.edgesN)){
                spN<-spnamesN[i] # get species names
                termabunN[i]<-(terminal.edgesN[spN])*(as.integer((ABUND_N_ONLY_S[spN][2,]))-1)
            }
            PAEdividendN<-(sum(unlist(termabunN)))+PDvalueN
            PAEdivisorN<-PDvalueN+(((sum(as.integer(ABUND_N_ONLY_S[2,][which(ABUND_N_ONLY_S[2,]>=1)]))/length(ABUND_N_ONLY_S[which(ABUND_N_ONLY_S[2,]>=1)]))-1)*sum(terminal.edgesN))
            PAEquotientN<-PAEdividendN/PAEdivisorN
            
            SIMULATION[rep,1] <- PAEquotient/PAEquotientN
            rep <- rep+1
            # p-value for SESPD is included in ses.pd output
        }
        p_v <- sum(SIMULATION > OVERALL_PAE_RATIO)/numReps
        if(p_v >.5){p_v <- 1-p_v}
        ################################################################################################################################
        ################################################################################################################################
        ################################################################################################################################
        # Plot the null dist
        d <- density(SIMULATION[,1])
        #to shade 95% area under curve
        q2_5 <- quantile(SIMULATION[,1], 0.025)
        q97_5 <- quantile(SIMULATION[,1], 0.975)
        x1 <- min(which(d$x > q2_5))
        x2 <- max(which(d$x < q97_5))
    
        #############################
        ###### EXPORT ###############
        #############################
        EXPORT_PAE <- rbind(EXPORT_PAE,cbind(as.character(CURRENT_SITE),
                                             length(spnames),length(spnamesN),
                                             PAEquotient_O,PAEquotientN_O,OVERALL_PAE_RATIO,
                                             PDvalue_0,PDvalueN_0,
                                             PAEdividend,PAEdividendN,
                                             PAEdivisor,PAEdivisorN))
        #EXPORT_SESPD <- rbind(EXPORT_SESPD,cbind(as.character(CURRENT_SITE),SESPD_Sim,SESPD_N_ONLY_Sim))
        EXPORT_SESPD_Groups <- rbind(EXPORT_SESPD_Groups,cbind(rbind(paste(CURRENT_SITE,"_NN",sep=""),paste(CURRENT_SITE,"_N",sep="")),SESPD_Groups))
    }
    # write each site 
    filename_PAE_S <- paste(DIRECTORY,as.character(CURRENT_SITE),"_PAE.csv",sep="")
    filename_SESPD_S <- paste(DIRECTORY,as.character(CURRENT_SITE),"_SESPD.csv",sep="")
    filename_SESPD_Groups_S <- paste(DIRECTORY,as.character(CURRENT_SITE),"_SESPD_Groups.csv",sep="")
    colnames(EXPORT_PAE) <- c("Site",
                                  "Num_Tips","Num_Tips_N_ONLY",
                                  "PAE_Quotient","PAE_Quotient_N_ONLY","PAE_Ratio",
                                  "PD_Value","PD_Value_N_ONLY",
                                  "PAE_Dividend","PAE_Dividend_N_ONLY",
                                  "PAE_Divisor","PAE_Divisor_N_ONLY")
    colnames(EXPORT_SESPD_Groups)[1] <- "Site"
    write.csv(EXPORT_PAE,filename_PAE_S,row.names=FALSE)
    write.csv(EXPORT_SESPD,filename_SESPD_S,row.names=FALSE)
    write.csv(EXPORT_SESPD_Groups,filename_SESPD_Groups_S,row.names=FALSE)
}
