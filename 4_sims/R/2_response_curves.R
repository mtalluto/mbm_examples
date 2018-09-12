library(mbm)

mbm1 <- readRDS("4_sims/res/mbm/mbmmod_25_1.rds")
gdm1 <- readRDS("4_sims/res/gdm_noperm/gdmmod_25_1.rds")

mbm2 <- readRDS("4_sims/res/mbm/svgpmod_225_1.rds")
gdm2 <- readRDS("4_sims/res/gdm_noperm/gdmmod_225_1.rds")

mbm3 <- readRDS("4_sims/res/mbm/svgpmod_2025_1.rds")
gdm3 <- readRDS("4_sims/res/gdm_noperm/gdmmod_2025_1.rds")


cex_pt <- 0.2
pch <- 16
yl = c(0, 1.1)
tcex = 0.9
cols <- c("#ff644b", "#005ab8", "#c260a5")
tp='png'

quartz(width=6.2, height=8.5, pointsize=10, type=tp, file=paste0("4_sims/sims.", tp), 
	bg='white', dpi=600)

par(mfrow=c(3,2), bty='n', mar=c(5,4,0.5,1))

rc(mbm1, ylim=yl, lwd=2, cex_pt = cex_pt, col_line=cols[1], xlab="Environmental Distance", 
	ylab="Sørensen (Taxonomic) Dissimilarity", pch=pch)
text(grconvertX(0.02, from='npc'), grconvertY(0.92, from='npc'), 'n = 25', pos=4, cex=tcex)
plot(gdm1$ecological, gdm1$observed, xlab = "GDM Ecological Distance", ylab = "", type = "n", 
	ylim=yl)
points(gdm1$ecological, gdm1$observed, pch = 16, cex = 0.2, col = "#666666")
overlayX <- seq(from = min(gdm1$ecological), to = max(gdm1$ecological), length = 200)
overlayY <- 1 - exp(-overlayX)
lines(overlayX, overlayY, lwd = 2, col=cols[2])


rc(mbm2, ylim=yl, lwd=2, cex_pt = cex_pt, col_line=cols[1], xlab="Environmental Distance", 
	ylab="Sørensen (Taxonomic) Dissimilarity", pch=pch)
text(grconvertX(0.02, from='npc'), grconvertY(0.92, from='npc'), 'n = 225', pos=4, cex=tcex)
plot(gdm2$ecological, gdm2$observed, xlab = "GDM Ecological Distance", ylab = "", type = "n", 
	ylim=yl)
points(gdm2$ecological, gdm2$observed, pch = 16, cex = 0.2, col = "#666666")
overlayX <- seq(from = min(gdm2$ecological), to = max(gdm2$ecological), length = 200)
overlayY <- 1 - exp(-overlayX)
lines(overlayX, overlayY, lwd = 2, col=cols[2])


rc(mbm3, ylim=yl, lwd=2, cex_pt = cex_pt, col_line=cols[1], xlab="Environmental Distance", 
	ylab="Sørensen (Taxonomic) Dissimilarity", pch=pch)
text(grconvertX(0.02, from='npc'), grconvertY(0.92, from='npc'), 'n = 2025', pos=4, cex=tcex)
plot(gdm3$ecological, gdm3$observed, xlab = "GDM Ecological Distance", ylab = "", type = "n", 
	ylim=yl)
points(gdm3$ecological, gdm3$observed, pch = 16, cex = 0.2, col = "#666666")
overlayX <- seq(from = min(gdm3$ecological), to = max(gdm3$ecological), length = 200)
overlayY <- 1 - exp(-overlayX)
lines(overlayX, overlayY, lwd = 2, col=cols[2])

dev.off()


