###################################################################
######## Invasive Binomial Builder ################################
###################################################################
Invasive.Binomial.Builder <- function(file){
    INVASIVE <- file
    IBinomial <- rep()
    for (name in INVASIVE$ScientificName){
        temp <- unlist(strsplit(name, split = "\\s+"))[1:2]
        name <- paste(temp, collapse="_")
        IBinomial <- rbind(IBinomial, name)
    }
    IBinomial <- data.frame(IBinomial)
    colnames(IBinomial) <- c("IBinomial")
    INVASIVE <- data.frame(INVASIVE$Symbol)
    colnames(INVASIVE) <- c("Invasive")
    EXPORT <- cbind(IBinomial,INVASIVE)
    return(EXPORT)
}
