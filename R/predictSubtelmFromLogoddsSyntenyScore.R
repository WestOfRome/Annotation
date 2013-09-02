y=read.table('chr4.tab', sep="\t")

nr=nrow(y)
fraction=30/nr
#print(fraction) 
#fraction=.1

start=20
stop=20
col=rainbow( (stop-start+1) )

par(new=FALSE)
ylim=c(-5,1)

for (i in start:stop) { 
	x=y[,i]
	q= na.omit( x[ abs(x) < Inf ] )

	low=lowess(q, f= fraction)
	max  = mean(low$y)
	low$y = log2( low$y / max )

	offset=i-start+1
	#plot(  low, col=col[offset], ylim=ylim )
	
#	plot( cumsum( c(head(y[,i],50),tail(y[,i],50)) ), col=col[offset] )
	plot( cumsum( y[,i]), col=col[offset] )
		
	mtext(i,line=-2*offset, col=col[offset], cex=2)
	par(new=TRUE)
}

abline(v=11)
abline(v=(nr - 18))

#plot( ecdf(z$V21 ) )
#plot(cumsum(y$V21))

