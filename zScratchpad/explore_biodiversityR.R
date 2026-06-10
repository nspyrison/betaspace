library(vegan)
library(MASS)
library(BiodiversityR)

data(dune)
distmatrix <- vegdist(dune)
Ordination.model1 <- NMSrandom(distmatrix, perm=100, k=2)
Ordination.model1 <- add.spec.scores(Ordination.model1, dune, 
                                     method='wa.scores')
Ordination.model1
