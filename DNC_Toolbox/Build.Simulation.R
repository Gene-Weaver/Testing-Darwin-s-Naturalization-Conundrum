# Generate Simulation Tree
# use the NONsite-specific 'tree'
Build.Simulation <- function(tree,SITE_INVASIVE_COUNT){
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
    return(list(ABUND_S,ABUND_N_ONLY_S,sim_tree_RAND,sim_tree_RAND_N_ONLY))
}