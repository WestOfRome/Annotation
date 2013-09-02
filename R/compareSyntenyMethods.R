#y=read.table('all3.tab', sep="\t", skip=6)
sub=y[,c(5:13)]

layout(matrix(c(1:4),nrow=2))

print(round(cor(sub, use='pairwise.complete.obs', meth='sp')*100))
image(round(cor(sub, use='pairwise.complete.obs', meth='sp')*100), col=rev(rainbow(10,start=0,end=.6)))

# plot( sub )
h.index=8
l.index=11
hyper=y[,h.index]
loss=y[,l.index]

print( cor.test(loss,hyper,meth="sp"));

plot(density(hyper), main="Hyper7_+YGOB")
abline(v=5,col="red")

plot(density(loss), main="LOSS3_-YGOB-SGD")
abline(v=40,col="red")

plot(loss,hyper)
abline(h=5,col="green")
abline(h=15,col="red")
abline(v=40,col="green")
abline(v=100,col="red")

loss.not.hyper=subset(y, y[,l.index] >=40 & y[,h.index] <= 5)
hyper.not.loss=subset(y, y[,l.index] <=40 & y[,h.index] >= 5)
loss.and.hyper=subset(y, y[,l.index] >=40 & y[,h.index] >= 5)

print(nrow(loss.and.hyper))
print(nrow(loss.not.hyper))
print(nrow(hyper.not.loss))

