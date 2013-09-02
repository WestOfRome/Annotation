################################################
# Copyright Devin Scannell 2013 
# Simple script to use YGOB ohnolog data to 
# develop sensible scoring and statistics for 
# the ohnolog2() algorithm.
################################################

options(scipen=3);
options(width=200);

color=c('red','blue');
ymax.def=0.15

file="ohno2_long-win-10.tab"
file="ohno2_delta-nu.tab"

################################################
# read in algorithm scoring data  
################################################

x <- read.table(file,header=T,skip=0, row.names=1)

# remove junk columns from the output 

x$X. <- NULL
x$X..1 <- NULL
backup<-x

# there are duplicate rows -- need to keep only the "best" one 

x<-backup
best.rows<-c()
for ( anc in unique(x$AncLocus) ) {
	x.anc <- subset(x, x$AncLocus==anc)
	best.score <- as.numeric(row.names(x.anc[order(x.anc$Score,decreasing=T),])[1])
	#print(c(anc,nrow(x.anc),best.score), quote=F, digits=1)
	
	for (rn in row.names(x.anc)) {
		if ( rn==best.score ) {
			best.rows<-c(best.rows,as.numeric(rn))
		}
	}
}

print(nrow(x))
y<-x[best.rows,]
print(nrow(y))
x<-y

################################################
# read in actual S. cerevisiae ohnologs from YGOB 
# and label rows in data to identify true positives 
################################################

ohno <- scan("ohnologs.ygob.scer.nrl",what="character");
m.x.ohno <- match(ohno, x$AncLocus)

x$TrueP <- 0
x$TrueP[m.x.ohno] <- 1
table(x$TrueP)

################################################
# look at outliers to find any issues ….  
################################################

options(width = 200)
sub <- subset(x, x$TrueP==0 & x$Score > 10)
nrow(sub)
print(sub)
sub <- subset(x, x$TrueP==1 & x$Score < -10)
nrow(sub)
print(sub)

################################################
# prep work before graphing 
################################################

align.min <- function (rx) {
	return(min(as.numeric(rx[14:15])))
}
x$Align_Max <- as.numeric(apply(x,1,align.min))
x$OhnoN<-x$OHNO/x$Align_Max
x$CrossX<-as.numeric(x$CROSS/x$SAME)

################################################
# look at data …
# Define function so we can recycle later 
################################################

plot.align.vars <- function (df, features, ymax, vline) {

	layout(matrix(1:length(features),nrow=2,byrow=T))
	for (feat in features) {
		
	sub<-subset(df,df[[feat]]<Inf & ! is.na(df[[feat]]))
	sp.coef<-cor.test(sub[[feat]],sub$TrueP)
	#print(c(sp.coef$estimate,sp.coef$p.value),digits=3);
		
	q1 <- na.omit(df[[feat]][df$TrueP==1])
	q0 <- na.omit(df[[feat]][df$TrueP==0])
	s1 <- q1[q1<Inf]
	s0 <- q0[q0<Inf]
	
	xl <- c(min(c(s1,s0)),max(c(s1,s0)))
	yl <- c(0,ifelse(feat=="OHNO", .5, ifelse(feat=="DELTA", .005, ymax) ))
	
	par(new=F)
	plot( 
		density(s0), 
		xlim=xl, ylim=yl, col="red", 
		xlab=feat, ylab="Density", main=feat 
	)
	par(new=T)
	plot( 
		density(s1), 
		xlim=xl, ylim=yl, col="blue", 
		xlab=feat, ylab="Density", main=feat
	)
	mtext(paste(round(sp.coef$estimate,3),round(sp.coef$p.value,12),sep=" / "), side=3, line=-2,col="blue")
	
	if (vline) {
		abline(v=vline,col="green")	
	}
	# numbres .. 

	form <- formula( paste(feat,"TrueP", sep=" ~ ") )
	print(aggregate(form, x, mean))
	}
}

features <- c('Score','SAME', 'CROSS', 'GAP','OHNO','DELTA') # 'CrossX', 'OhnoN',
plot.align.vars(x, features, ymax.def)

file.eps=paste(file,'eps',sep='.')
dev.copy2eps(file=file.eps);

################################################
# look at data …
################################################

var='CROSS'
rx=c(1,3)
options(width = 200)
sub <- subset(x, x$TrueP==0 & x[[var]] > rx[3])
nrow(sub)
print(sub)
sub <- subset(x, x$TrueP==1 & x[[var]] < rx[1])
nrow(sub)
print(sub)

################################################
# Plot OHNO and DELTA 
################################################

var1='OHNO'
var2='CROSS'
#layout(matrix(1:2,nrow=1))
layout(1)
xl=range(x[[var1]])
yl=range(x[[var2]])

for (i in c(0,1)) {
	par(new=ifelse(i==0,FALSE,TRUE))
	#print(nrow(subset(x,x$TrueP==i)),col=color[i+1])

	plot( 
	x[x$TrueP==i,var1], x[x$TrueP==i,var2],
	col=color[i+1], xlim=xl, ylim=yl,
	xlab=var1, ylab=var2, main=i
	)
}

cx <- x[x$TrueP==i,c('OHNO','DELTA')]

################################################
# make score 
################################################

ymax=.1
cutoff.score=2

# tests <- c('test1','test2', 'test3', 'test4')
#x$test1 <- x$Score
#x$test3 <- x$test2 + x$KC
#x$test4 <- x$test3 - 2*x$OHNO + 2*x$CROSS  

test <- c('Score','test', 'OHNO', 'CROSS')
x$test <- x$OHNO + x$CROSS

print(summary(x$Score))
print(summary(x$KC))

plot.align.vars(x,test.x,ymax, cutoff.score)

x$BestGuess <- ifelse(x$test2 >= cutoff.score, "Ohno", "No Oh")

for (feat in tests) {
	x[[feat]]<-NULL
}

print(file);
table(x$TrueP,x$BestGuess)

file.pdf=paste(file,'pdf',sep='.')
if ( file.exists(file.pdf) ) {
	file.remove(file.pdf)
}
dev.copy2pdf(file=file.pdf);
