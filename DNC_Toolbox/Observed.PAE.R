# Calculate Observed PAE all and N-only
# VALUES TO EXPORT MARKED WITH ***^***
#
# Observed PAE -> ses.pd gets 1000 iterations but only 10 for the simulations
# PD is phylogenetic distance by summing edge lengths, SESPD is Phylogenetic diversity
# SESPD_Groups calculates sespd values based on native, NN status i.e. for 73 N an then 10 NN (Not very helpful for our purposes)
#   so it's calculated once per tree per site

# PAEdividend = numerator of PAE
# PAEdivisor = denominator of PAE
# PAEquotient = PAE
Observed.PAE <- function(TREE,TREE_N_ONLY,COMM,ABUND,ABUND_N_ONLY){
    print("Observed PAE")
    terms <- TREE$edge[,2] <= Ntip(TREE)
    terminal.edges <- TREE$edge.length[terms]
    names(terminal.edges)<-TREE$tip.label
    
    # Sum branch lengths for phylo dist
    PDvalue<-sum(TREE$edge.length)
    PDvalue_0 <- PDvalue # ***^***
    
    # SESPD All Taxa
    SESPD_Groups <- ses.pd(COMM,TREE,null.model = "taxa.labels", runs = 999, iterations = 1000) # ***^***
    SESPD <- ses.pd(ABUND,TREE,null.model = "taxa.labels", runs = 999, iterations = 1000)
    # This SESPD value is for the whole site tree, will be compared to the SESPD_N_ONLY values
    SESPD <- SESPD[1,] # ***^***
    
    #Calculate PAE
    spnames<-names(terminal.edges)
    termabun<-numeric(length(terminal.edges))
    for(i in 1:length(terminal.edges)){
        sp<-spnames[i] # get species names
        termabun[i]<-(terminal.edges[sp])*(as.integer((ABUND[sp][2,]))-1)
    }
    PAEdividend<-(sum(unlist(termabun)))+PDvalue_0
    PAEdivisor<-PDvalue_0+(((sum(as.integer(ABUND[2,][which(ABUND[2,]>=1)]))/length(ABUND[which(ABUND[2,]>=1)]))-1)*sum(terminal.edges))
    PAEquotient<-PAEdividend/PAEdivisor
    #print(PAEquotient)
    
    #######################################
    ####### Calculate PAE - N Only ########
    #######################################
    termsN <- TREE_N_ONLY$edge[,2] <= Ntip(TREE_N_ONLY)
    terminal.edgesN <- TREE_N_ONLY$edge.length[termsN]
    names(terminal.edgesN)<-TREE_N_ONLY$tip.label
    PDvalueN<-sum(TREE_N_ONLY$edge.length)
    PDvalueN_0 <- PDvalueN # ***^***
    
    SESPD_N_ONLY <- ses.pd(ABUND_N_ONLY,TREE_N_ONLY,null.model = "taxa.labels", runs = 999, iterations = 1000)
    SESPD_N_ONLY <- SESPD_N_ONLY[1,] # ***^***
    
    spnamesN<-names(terminal.edgesN)
    termabunN<-numeric(length(terminal.edgesN))
    for(i in 1:length(terminal.edgesN)){
        spN<-spnamesN[i] # get species names
        termabunN[i]<-(terminal.edgesN[spN])*(as.integer((ABUND_N_ONLY[spN][2,]))-1)
    }
    PAEdividendN<-(sum(unlist(termabunN)))+PDvalueN_0
    PAEdivisorN<-PDvalueN_0+(((sum(as.integer(ABUND_N_ONLY[2,][which(ABUND_N_ONLY[2,]>=1)]))/length(ABUND_N_ONLY[which(ABUND_N_ONLY[2,]>=1)]))-1)*sum(terminal.edgesN))
    PAEquotientN<-PAEdividendN/PAEdivisorN
    #print(PAEquotientN)
    
    OVERALL_PAE_RATIO <- PAEquotient/PAEquotientN
    
    EXPORT_SESPD <- cbind(SESPD,SESPD_N_ONLY) #1
    EXPORT_SESPD_GROUPS <- SESPD_Groups #2
    EXPORT_PAE_ALL <- cbind(PDvalue_0,PAEdividend,PAEdivisor,PAEquotient,OVERALL_PAE_RATIO) #3
    EXPORT_PAE_N_ONLY <- cbind(PDvalueN_0,PAEdividendN,PAEdivisorN,PAEquotientN) #4
    return(list(EXPORT_SESPD,EXPORT_SESPD_GROUPS,EXPORT_PAE_ALL,EXPORT_PAE_N_ONLY,OVERALL_PAE_RATIO))
}





