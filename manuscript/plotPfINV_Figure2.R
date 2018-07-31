#GeneVsGene Plot
geneTop <- 'PF3D7_0402500'
geneBot <- 'PF3D7_0402700'

bpTop <- 155877:156984 # Base positions for top gene
bpBot <- 162511:163724 # Base positions for bottom gene

exonsTop <- c(155878:155947,156068:156983)
exonsBot <- c(162512:163016,163020:163401,163405:163438,163657:163723)

intronsTop <- bpTop[!(bpTop %in% exonsTop)]
intronsBot <- bpBot[!(bpBot %in% exonsBot)]

Top2Bot <- (bpBot[1]-bpTop[1])

topINV <- c(356, 1108) # Matching region in top gene
botINV <- c(1,757) # Matching region in bottom gene

TopAsBot <- bpTop + (bpBot[1]-bpTop[1])
spacing <- 200
diff <- 4 # Calculate this, is displacement of bp5 and bp7
topticksPos <- TopAsBot[TopAsBot %% spacing == 0] + diff
topticksLab <- bpTop[TopAsBot %% spacing == 0] + diff

topVar <- list(c(156104,156665), c(156265,156665))
botVar <- list(c(162832,163401), c(162832,163234))
colINVs <- c('red','blue')

h <- 4
w <- 8
u <- 'in'
r <- 300
tiff('PfINV_Figure2.tiff', height=h, width=w, unit=u, res=r)
par(cex=0.7)
#par(mfrow=c(3,1))

# Plot base genes
plot(bpBot, rep(1,length(bpBot)), col='white', yaxt='n', ylab='', xlab='', bty='n', ylim=c(0.75,1.25)) #base
points(bpBot[1:length(bpTop)], rep(1.2,length(bpTop)), type='l', lwd=1) #5
points(bpBot, rep(0.8,length(bpBot)), type='l', lwd=1) #7
points(exonsTop+Top2Bot, rep(1.2,length(exonsTop)), type='p', lwd=2, pch=15) #5
points(exonsBot, rep(0.8,length(exonsBot)), type='p', lwd=2, pch=15) #7
axis(side=3,at=topticksPos, labels=topticksLab)
text(min(bpBot), c(1.22, 0.82), labels=c(geneTop, geneBot), pos=4)
legend('left', legend=c('BLAST (81.6%)','INV018-INV021','INV019-INV020'), col=c('black', colINVs), lty=c(NA,1,1), pch=c(22,NA,NA))
legend('left', legend=c('','',''), col=c(rgb(0,1,0,0.2), colINVs), lty=c(NA,1,1), pch=c(15,NA,NA), bty='n')

# Plot gene-gene INV coordinates
polygon(c(bpBot[botINV[1]],bpBot[botINV[2]],bpBot[topINV[1]], bpBot[topINV[2]]), c(1.2,1.2,0.8,0.8), col=rgb(0,1,0,0.2))

# Plot INV variants positions
#polygon(c(topVar[[1]]+(bpBot[1]-bpTop[1]),botVar[[1]]), c(1.2,1.2,0.8,0.8), lty=2, col=rgb(0,0,1,0.2))
for (x in 1:length(topVar)) {
	# plot(bpBot, rep(1,length(bpBot)), col='white', yaxt='n', ylab='', xlab='', bty='n', ylim=c(0.75,1.25)) #base
	# points(bpBot[1:length(bpTop)], rep(1.2,length(bpTop)), type='l', lwd=2) #5
	# points(bpBot, rep(0.8,length(bpBot)), type='l', lwd=2) #7
	# axis(side=3,at=topticksPos, labels=topticksLab)
	# text(min(bp7), c(1.215, 0.815), labels=c(geneTop, geneBot), pos=4)
	# polygon(c(topVar[[x]]+(bpBot[1]-bpTop[1]),botVar[[x]]), c(1.2,1.2,0.8,0.8), lty=2, col=rgb(0,0,1,0.2))
	lines(topVar[[x]]+Top2Bot, rep(1.215+x/100, 2), col=colINVs[x], lwd=2)
	lines(botVar[[x]], rep(0.785-x/100, 2), col=colINVs[x], lwd=2)
}

dev.off()
