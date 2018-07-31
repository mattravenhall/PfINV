library(RColorBrewer)

fig_label <- function(text, region="figure", pos="topleft", cex=1.7, ...) {
 
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

dat <- read.csv('Pacbio_INV_wInHouse.csv')
#dat <- read.csv('Pacbio_INV_clean.csv')

# Size
h <- 8
w <- h*1.5
u <- 'in'
r <- 120

xdivide <- 0.4
xdivideB <- 0.7

ydivide <- 0.55
ydivideA <- 0.3
ydivideB <- 0.6
#ydivideC <- 0

chrLengths <- c(643292,947102,1060087,1204112,1343552,1418244,1501717,1439563,1541723,1687655,2038337,2271478,2895605,3291871)

defaultMAI <- c(1.02,0.82,0.82,0.42)

tiff('PfINV_Figure1.tiff', height=h, width=w, unit=u, res=r)
par(cex=0.7)

# Plot chromosome locations #######
par(fig=c(0,xdivide,0,1), mai=c(0.7,0.7,0.3,0.1))
bp <- barplot(chrLengths, names=1:length(chrLengths), ylab='Chromosome', xlab='Position (bp)', xaxt='n', horiz=T, las=1, col='lightgrey', border='white')
plotINVs <- function(df) {
  chr <- as.numeric(df[which(names(df) == 'Chr')])
  start <- df[which(names(df) == 'Start_3D7')]
  end <- df[which(names(df) == 'End_3D7')]
  rect(start, bp[chr,]-0.5, end, bp[chr,]+0.5, col='black', lwd=0.5, lend=1, border='black')
}
axis(1, at=seq(0,3000000, 1000000), labels = formatC(seq(0,3000000, 1000000), big.mark = ",", format = "d"), las = 1)
apply(dat, 1, plotINVs)
fig_label('A.')
box("figure")
###################################

# Plot INV type pie ###############
par(fig=c(xdivide,1,ydivide,1), new=TRUE, mai=c(0.1,0.1,0.2,0.1)) #x1, x2, y1, y2
INVtypes <- sort(table(dat$Type),ascending=FALSE)
pieCols <- rev(brewer.pal(length(INVtypes)+1, "Greens")) #rev(brewer.pal(length(INVtypes), "Accent")) #rainbow(length(INVtypes), s=0.5)
pie(INVtypes, labels=NA, clockwise=TRUE, col=pieCols, radius=1, border='white')
#legend(0.2,-0.5, rev(names(table(dat$Type))), cex=1, fill=pieCols)
legend('right', rev(names(table(dat$Type))), cex=1, fill=pieCols)
fig_label('B.')
box("figure")
###################################

## Plot length group counts #######
par(fig=c(xdivide,xdivideB,0,ydivide), new=TRUE, mai=c(0.7,0.7,0.3,0.1))
dividers <- c(0,500,1000,5000,10000,999999999999)
groupNames <- c('0-500 bp','500-1,000 bp','1,000-5,000 bp','5,000-10,000 bp', '>10,000 bp')
lengthCols <- rev(brewer.pal(length(groupNames)+1, "Reds"))
values <- cut(dat$Len_3D7,dividers)
heights <- as.vector(table(cut(dat$Len_3D7,c(0,500,1000,5000,10000,999999999999))))
bars <- plot(values, names=groupNames, ylim=c(0,120), col=lengthCols, xlab='Length', ylab='Count', border='white')
text(bars, heights, paste('n = ', heights, sep=''), pos=3)
fig_label('C.')
box("figure")
###################################

# Plot identity groups ############
par(fig=c(xdivideB,1,0,ydivide), new=TRUE, mai=c(0.7,0.7,0.4,0.1))
#identityCols <- c('brown', 'grey', 'goldenrod')
#barplot(table(cut(dat$Identity,c(0,90,95,100))), col=identityCols, names=c('<90%','90-95%','>95%'), xlab='Identity', ylab='Count', ylim=c(0,80))
dividers <- c(0,90,95,100)
groupNames <- c('<90%','90-95%','>95%')
identityCols <- rev(brewer.pal(length(groupNames)+1, "Blues")) #c('brown', 'grey', 'goldenrod')
values <- cut(dat$Identity,dividers)
heights <- as.vector(table(cut(dat$Identity,dividers)))
bars <- plot(values, names=groupNames, ylim=c(0,100), col=identityCols, xlab='Identity', ylab='Count', border='white')
text(bars, heights, paste('n = ', heights, sep=''), pos=3)
fig_label('D.')
box("figure")
###################################
dev.off()


#############################################################################

if (FALSE) {
library(beanplot)
# Plot intergenic & highly variable INV lengths vs other
tiff('InterAndVariableLengthsVsOthers.tiff', height=h, width=h*2, res=r, unit=u)
x <- dat$Len_3D7[dat$Type %in% c('Intergenic', 'Highly variable region')]
y <- dat$Len_3D7[!(dat$Type %in% c('Intergenic', 'Highly variable region'))]
beanplot(x,y, col=c('lightgrey','black'), horizontal=TRUE, names=c('Intergenic &\nHighly Variable', 'Other'), 
	xlab='Length (bp)', what=c(F,T,F,T), ll=0.05, main='Comparison of Inversion Lengths')
dev.off()
}

#############################################################################
