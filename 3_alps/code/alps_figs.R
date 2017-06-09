library(mbm)
taxModel <- readRDS("3_alps/res/taxModel.rds")
taxModelMF <- readRDS("3_alps/res/taxModelMF.rds")
funModel <- readRDS("3_alps/res/funModel.rds")
phyModel <- readRDS("3_alps/res/phyModel.rds")

## FIG 2 - look at response curves of models
tp <- "pdf"
cex_pt <- 0.2
pch <- 16

quartz(width=6.5, height=2.7, pointsize=11, type=tp, file=paste0("3_alps/img/rc.", tp), bg='white', dpi=600)
par(mfrow=c(1,3), bty='n', mar=c(5,4,0.5,1), bty='n')
cols <- c("#ff644b", "#005ab8")
rc(taxModel, ylim=c(0,1), lwd=2, cex_pt = cex_pt, col_line=cols[1], xlab="", ylab="Sørensen (Taxonomic) Dissimilarity", pch=pch)
rc(taxModelMF, add=TRUE, col_pt = NA, col_line = cols[2], lwd=2)
legend('bottomleft', lwd=2, legend=c("Mean function = 0", "Increasing mean function"), col=cols, bty='n')

rc(funModel, lwd=2, cex_pt = cex_pt, col_line=cols[1], xlab="", ylab="Functional MPD", pch=pch)
rc(phyModel, lwd=2, cex_pt = cex_pt, col_line=cols[1], xlab="", ylab="Phylogenetic MPD", pch=pch)
mtext(side=1, "Environmental Distance", outer=TRUE, line=-1, cex=0.9)
dev.off()


