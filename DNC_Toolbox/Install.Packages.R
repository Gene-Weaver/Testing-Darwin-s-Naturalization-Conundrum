Install.Packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("plyr",
              "dplyr",
              "ape",
              "phytools",
              "tcltk",
              "picante",
              "evobiR",
              "phyloch",
              "tictoc")
update.packages(ask = FALSE, checkBuilt = TRUE)
Install.Packages(packages)

