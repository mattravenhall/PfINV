fig_label <- function(text, region="figure", pos="topleft", cex=1.2, ...) {
 
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
 
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
 
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
 
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
 
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
 
  x1 <- switch(pos,
    topleft     =x[1] + sw, 
    left        =x[1] + sw,
    bottomleft  =x[1] + sw,
    top         =(x[1] + x[2])/2,
    center      =(x[1] + x[2])/2,
    bottom      =(x[1] + x[2])/2,
    topright    =x[2] - sw,
    right       =x[2] - sw,
    bottomright =x[2] - sw)
 
  y1 <- switch(pos,
    topleft     =y[2] - sh,
    top         =y[2] - sh,
    topright    =y[2] - sh,
    left        =(y[1] + y[2])/2,
    center      =(y[1] + y[2])/2,
    right       =(y[1] + y[2])/2,
    bottomleft  =y[1] + sh,
    bottom      =y[1] + sh,
    bottomright =y[1] + sh)
 
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
 
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

plotAlignment <- function(givenRegion,regionName, alignA, alignB, alignC, linestyle=c(1,1), fig='A.', boxCol='goldenrod') {
  # Desired region
#  chr <- 12
#  bpA <- 968000
#  bpB <- 980000
#  regionName <- 'gch1'
  
  # Alignments to plot
#  alignA <- c(967400,976400)
#  alignB <- c(976400,971100)
#  alignC <- c(971100,982000)

  chr <- givenRegion[1]
  bpA <- givenRegion[2]
  bpB <- givenRegion[3]
  
  # Annotation
  annotation <- read.csv('Pf3D7_annotation.tsv', sep='\t')
  annot <- annotation[(annotation$Chromosome == chr) & (annotation$Start <= bpB) & (annotation$End >= bpA),]
  
  # Setup
  region <- c(bpA,bpB)
  alignA[1] <- region[1]
  alignC[2] <- region[2]
  alignments <- c(alignA, alignB, alignC)
  spacer <- (region[2] - region[1]) / 100
  
  ## 
  baseheights <- rep(max(region),1+length(alignments)/2)
  par(mar=c(2,2,2,2))
  bp <- barplot(baseheights, horiz=T, col='white', border=NA, xlab='', main=regionName, xlim=region, cex.axis=0.7)
  
  plotGene <- function(x) {
  	cols <- names(x)
  	start <- as.numeric(x[which(cols == 'Start')])
  	end <- as.numeric(x[which(cols == 'End')])
  	gene <- x[which(cols == 'Feature')]
  	rect(start, bp[1,]-0.5, end, bp[1,]+0.5, col='grey', border=NA)
  	text((start+end)/2, bp[1,], labels=gene, cex=0.75)
  }
  apply(annot, 1, plotGene)
#  lines(region, rep(bp[1,]+0.5,2), lty=2)
  
  ## Plot alignment tracks
  rect(alignA[1], bp[4,]-0.5, alignA[2], bp[4,]+0.5, col=boxCol, border=NA)	# start
  text((alignA[1]+alignA[2])/2, bp[4,], labels='>')
  rect(alignB[1], bp[3,]-0.5, alignB[2], bp[3,]+0.5, col=boxCol, border=NA)	# join
  text((alignB[1]+alignB[2])/2, bp[3,], labels='<')
  rect(alignC[1], bp[2,]-0.5, alignC[2], bp[2,]+0.5, col=boxCol, border=NA)	# end
  text((alignC[1]+alignC[2])/2, bp[2,], labels='>')
  
  # # Plot join lines
  # # Before first
  lines(c(alignA[1]-(10*spacer), alignA[1]), rep(bp[4,],2))
  # # For each gap between
  lines(c(alignA[2], max(alignA[2],alignB[1])+spacer), rep(bp[4,],2), lty=linestyle[1]); lines(rep(max(alignA[2],alignB[1])+spacer,2), c(bp[4,], bp[3,]), lty=linestyle[1]); lines(c(max(alignA[2],alignB[1])+spacer,alignB[1]), rep(bp[3,],2), lty=linestyle[1])
  lines(c(min(alignB[2],alignC[1])-spacer, alignB[2]), rep(bp[3,],2), lty=linestyle[2]); lines(rep(min(alignB[2],alignC[1])-spacer,2), c(bp[3,], bp[2,]), lty=linestyle[2]); lines(c(min(alignB[2],alignC[1])-spacer,alignC[1]), rep(bp[2,],2), lty=linestyle[2])
  # # After last
  lines(c(alignC[2], alignC[2]+(10*spacer)), rep(bp[2,],2))
  fig_label(fig)
  box("figure")
}

# Set up outfile
h <- 5
w <- 8
u <- 'in'
r <- 300
tiff('PfINV_Figure3.tiff', height=h, width=w, unit=u, res=r)
options(scipen=999999)
par(mfrow=c(2,1))

plotAlignment(c(12,970000,978000),'gch1', c(967400,976400), c(976400,971100), c(971100,982000), fig='A.', boxCol='goldenrod')
plotAlignment(c(5,408000,420000),'pi4k', c(0,417648), c(417648,410463), c(410268,999999), linestyle=c(1,2), fig='B.', boxCol='maroon')

dev.off()