# Simulated MNTD

Simulated.MNTD <- function(TREE,ABUND,INVASIVE_DATA){
    temp_out <- Status.0.1.Builder(ABUND,INVASIVE_DATA)
    COMM <- temp_out[[1]]
    
    comm <- COMM
    phydist <- cophenetic(TREE)
    
    # MNTD and MNTD_MOD
    mntd_reg <- mntd(comm,phydist)
    mntd_reg <- mntd_reg[1]
    mntd_mod <- MNTD_MOD(comm, phydist, abundance.weighted = FALSE)
    mntd_mod <- mntd_mod[1]
    mntd_mod_ratio <- mntd_mod/mntd_reg
    
    MNTD_temp <- mntd(comm, phydist)    
    MNTD_temp[2]->MNTDnative
    MNTDbetween <- comdist(comm, phydist)
    MNTDratio <- MNTDbetween/MNTDnative
    
    #SES_MNTD_MOD
    SES <- SES_MNTD_MOD(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
    SES_reg <- ses.mntd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
    SES_MNTD_MOD_ratio <- SES[1,2]/SES_reg[1,2]
    
    #there appears to be a sign error in the ses_mntd_mod function... hence the lack of -1*sd RESOLVED=> we have the -1 , they don't
    SES_MOD_sd <- -1*((SES[1,2]-SES[1,3])/SES[1,4])
    SES_REG_sd <- -1*((SES_reg[1,2]-SES_reg[1,3])/SES_reg[1,4])
    SES_sd_ratio <- SES_MOD_sd/SES_REG_sd
    
    #OUTPUT <- cbind(rbind(cbind(mntd_reg,mntd_mod,MNTDratio,mntd_mod_ratio,SES_MNTD_MOD_ratio,SES_reg[1,2],SES[1,2],SES_MOD_sd,SES_REG_sd,SES_sd_ratio),cbind(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)),SES,SES_reg)
    #colnames(ALL_OUT) <- c("MNTD reg. (A)","MNTD MOD (B)","MNTD ratio (C)","MNTD MOD ratio (D)","SES MNTD MOD ratio (E)","SES MNTD reg (I)","SES MNTD MOD (J)","SES MOD sd (G)","SES reg sd (H)","SES_sd_ratio [SES_MOD_sd/SES_REG_sd] (F)","START SES MOD n_taxa","mntd.obs","mntd.rand.mean","mntd.rand.sd","mntd.obs.rank","mntd.obs.z","mntd.obs.p","runs",
    #                       "START SES reg n_taxa","mntd.obs","mntd.rand.mean","mntd.rand.sd","mntd.obs.rank","mntd.obs.z","mntd.obs.p","runs")
    
    
    return(SES_MNTD_MOD_ratio)
}



