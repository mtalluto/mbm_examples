taxBeta <- read.csv('3_alps/dat/tax-beta.csv')

taxMod <- mbm(taxBeta[,4:8], taxBeta$sor)