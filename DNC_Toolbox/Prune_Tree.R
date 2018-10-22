###
### Prune tree tips based on an imported .csv file
###

Prune.Tree <- function(SPECIES_TO_REMOVE,TREE,SAVE_FILE){
    PRUNED_TREE<-drop.tip(TREE,TREE$tip.label[SPECIES_TO_REMOVE])
    write.tree(PRUNED_TREE,file=SAVE_FILE)
    return(PRUNED_TREE)
}
