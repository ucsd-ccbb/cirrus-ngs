#install.packages("/shared/workspace/software/R-packages/DESeq", repos = NULL, type="source")
#install.packages("/shared/workspace/software/R-packages/limma", repos = NULL, type="source")
#install.packages("/shared/workspace/software/R-packages/qvalue", repos = NULL, type="source")

 .libPaths();

library(DESeq)
library(limma)
library(qvalue)

myfold <- function(x) {
   posx<-x>=0
   fc<-x
   fc[posx]<-2^x[posx]
   fc[!posx]<-(-1)*2^(-x[!posx])
   return(fc)
}


#infile<-commandArgs(trailingOnly = T)
args<-commandArgs(trailingOnly = T)
infile<-args[1]
rootpath<-args[2]
META<-read.table(infile,sep="\t",header=F) #read metadata
filename<-as.character(META[1,1])
ntext<-as.numeric(META[1,2])
pink<-0.01
red<-0.0001
mmx<-5 #y range of MA plots

k<-nrow(META)-1

print(paste("Doing a ",k,"-way analysis",sep=""))

conditions<-t(as.character(META[2:(k+1),1]))
nc<-length(conditions)
npairs<-nc*(nc-1)/2
outfile<-paste(c(rootpath,c(conditions,"lfdr.txt")),sep="",collapse="")

nrep<-t(META[2:(k+1),2])

nrto<-nrep
for (i in 2:k) {
   nrto[i]<-nrto[i-1]+nrep[i]
}

X<-read.table(filename,sep="\t",header=T) #read counts
nsamples<-sum(nrep) #number of samples

#keep only detected genes
rowsum<-apply(X[,(ntext+1):(ntext+nsamples)],1, sum)
detected<-rowsum>=(nsamples*16) #an arbitrary cutoff for sanity
print("Arbitrary mean count cutoff for sanity set to 16")
ngenes<-sum(detected)

names<-X[detected,1:ntext]
S<-round(X[detected,(ntext+1):(ntext+nsamples)])


group<-factor(c(rep(conditions,nrep)),levels=conditions)
#design<-model.matrix(~group) #for voom type 1
design<-model.matrix(~0+group) #for voom type 2
colnames(design)<-conditions

#----------------begin DESeq normalization---------------
expDesign<- data.frame(row.names=colnames(S), condition = c(rep(conditions,nrep)))
cds<-newCountDataSet(S, expDesign)
cds<- estimateSizeFactors(cds)
S_norm<-counts(cds, normalized=TRUE) #this is the "normalized" counts for plotting maybe
print("Relative library sizes:")
sizeFactors(cds)
#-----------------DESeq normalization done-----------------


S_norm<-pmax(S_norm,0)
data<-log2(S_norm+1/2)
fit<-lmFit(data,design)

contrast.matrix<-matrix(0,nrow=nc,ncol=npairs)
rownames(contrast.matrix)<-conditions
colnames(contrast.matrix)<-rep("",npairs)
l<-1
for (i in 1:(nc-1)) {
   for (j in (i+1):nc) {
      #doing the i,j pair of conditions
      contrast.matrix[i,l]<- -1
      contrast.matrix[j,l]<- 1
      colnames(contrast.matrix)[l]<-paste(conditions[j],"-",conditions[i],sep="")
      l<-l+1 #next column
   }
}

fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
p<-fit2$F.p.value #ANOVA p-value, needs to be adjusted
#qobj<-qvalue(p,pi0.meth="bootstrap")
#q<-qobj$qvalues #q-values 
# p will be converted into local fdr using Storey's method
lfdr0<-lfdr(p,pi0.method="bootstrap",transf="logit") #het lfdr from k-way p-values

locfdr<-1/(1+exp(fit2$lods)) #locfdr is a matrix with as many columns as contrasts
difstring<-colnames(contrast.matrix) #has a -
ratstring<-gsub("-","/",difstring) #has a /
#ratstring<-gsub("\\+","*",ratstring) #has a /
colnames(locfdr)<-paste("lfdr ",difstring,sep="")

log2r<-fit2$coefficients
colnames(log2r)<-paste("log2 ",ratstring,sep="")

fc<-myfold(log2r)
colnames(fc)<-paste("fc ",ratstring,sep="")

levels<-matrix(nrow=ngenes, ncol=k)
if (nrto[1] > 1) {
   levels[,1]<-apply(S_norm[,1:nrto[1]], 1, mean)
} else {
   levels[,1]<-S_norm[,1]
}
for (i in 2:k) {
   levels[,i]<-apply(S_norm[,(nrto[i-1]+1):nrto[i]], 1, mean)
}
colnames(levels)<-conditions

S_norm<-round(S_norm)
levels<-round(levels)
res0<-cbind(names,S_norm,levels,log2r,fc,locfdr)
if (k==2) {
o<-order(locfdr[,1])
} else {
o<-order(lfdr0)
res0<-cbind(res0,p,lfdr0)
}

write.table(res0[o, ], outfile, sep="\t", row.names=FALSE,col.names=TRUE,quote=FALSE)

#----------------------------------------------------------
pdffile<-paste(c(rootpath,c(conditions,"pca.pdf")),sep="",collapse="")
pdf(pdffile,paper="letter")
#par(pty="s")
par(xpd=T, mar=par()$mar+c(0,0,0,6), pty="s")
#col0<-c("red","green","blue")
col0<-c("red","green","blue","black","magenta","cyan","pink","lightgreen","lightblue","grey","orange","yellow")
color<-rep(col0[1:k],times=nrep)
txt<-as.character(1:nsamples)

pca.res <- prcomp(t(data))

xy<-pca.res$x[,1:2]
rge<-range(xy)
plot(xy,col=color,pch=16,cex=2,xlim=rge,ylim=rge,main="principal component analysis",xlab="PC 1",ylab="PC 2")
for (i in 1:nsamples) {
points(xy[i,1],xy[i,2],col=color[i],pch=16,cex=2)
text(xy[i,1],xy[i,2],txt[i],cex=.6,col="white")
}
legend(max(rge)*1.15,max(rge),conditions,pch=16,cex=1,col=col0[1:k])

w<-apply(data, 1, var)
d<-matrix(nrow = nsamples, ncol = nsamples)
for (i in 1:nsamples) {
for (j in 1:nsamples) {
d[i,j]<-sum((data[,i]-data[,j])^2*w) }}
xy<-cmdscale(d,k=2,eig=F)
plot(xy,pch=16,cex=2,col=color,main="principal coordinate analysis",xlab="distance",ylab="distance",xlim=range(xy),ylim=range(xy))
for (i in 1:nsamples) {
points(xy[i,1],xy[i,2],col=color[i],pch=16,cex=2)
text(xy[i,1],xy[i,2],txt[i],cex=.6,col="white")
}
legend(max(xy)*1.15,max(xy),conditions,pch=16,cex=1,col=col0[1:k])
dev.off()


#----------------------------------------------------------
pdffile<-paste(c(rootpath,c(conditions,"XY_scatterplots.pdf")),sep="",collapse="")
pdf(pdffile,paper="letter")
par(pty="s",mfrow=c(2,3))

labels<-conditions

#r<-matrix(runif(ngenes*nsamples), ngenes, nsamples)
e<-log2(levels)
q<-locfdr[,1]

for (i in 1:(k-1)) {
   for (j in (i+1):k) {
      x<-as.numeric(e[,i])
      y<-as.numeric(e[,j])
      plot(x,y,pch=16,cex=0.15,xlim=c(0,20),ylim=c(0,20),xlab=labels[i],ylab=labels[j])
      lines(c(0,20),c(0,20),col="blue")

      x_pink<-x[q<pink]
      y_pink<-y[q<pink]
      points(x_pink,y_pink,pch=16,cex=0.2,col="pink")
      x_red<-x[q<red]
      y_red<-y[q<red]
      points(x_red,y_red,pch=16,cex=0.2,col="red")
   }
}

#----------------------------------------------------------
dev.off()
#----------------------------------------------------------
pdffile<-paste(c(rootpath,c(conditions,"MA_scatterplots.pdf")),sep="",collapse="")
pdf(pdffile,paper="letter")
par(pty="s",mfrow=c(2,3))

labels<-conditions

#r<-matrix(runif(ngenes*nsamples), ngenes, nsamples)
e<-log2(levels)

for (i in 1:(k-1)) {
   for (j in (i+1):k) {
      x<-as.numeric(e[,i])
      y<-as.numeric(e[,j])
      a<-(x+y)/2
      m<-y-x
      m[m>mmx]<-mmx
      m[m<(-mmx)]<-(-mmx)
      plot(a,m,pch=16,cex=0.15,xlim=c(0,20),ylim=c(-mmx,mmx),xlab="log2 read count",ylab=paste("log2(",labels[j],"/",labels[i],")",sep="",collapse=""))
      lines(c(0,20),c(0,0),col="blue")

      a_pink<-a[q<pink]
      m_pink<-m[q<pink]
      points(a_pink,m_pink,pch=16,cex=0.2,col="pink")
      a_red<-a[q<red]
      m_red<-m[q<red]
      points(a_red,m_red,pch=16,cex=0.2,col="red")
   }
}
#----------------------------------------------------------
dev.off()

