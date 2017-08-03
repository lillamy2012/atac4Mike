args = commandArgs(trailingOnly=TRUE)


## first arg is type

type=args[1]

## rest files

files = args[2:length(args)]

### load libraries
library(wesanderson)
library(Cairo)

colors=c("red","blue")

is=list()
for (i in files){
  is[[i]] = read.table(i)
  is[[i]]$V1 = is[[i]]$V1-8 ## 
  rownames(is[[i]])=is[[i]]$V1
  is[[i]]=is[[i]][,2,drop=FALSE]
}
print(str(is))
## to matrix
nrows = max(sapply(is,function(x)  max(as.numeric(rownames(x)))))
ncols = length(is)
insertCounts = matrix(0,ncol=ncols,nrow=nrows)
colnames(insertCounts)=names(is)
rownames(insertCounts)=1:nrows
print(head(insertCounts))

for (i in names(is)){
  insertCounts[rownames(is[[i]]),i]=unlist(is[[i]])
}

print(t(t(insertCounts[1:250,])/colSums(insertCounts)))

#create plots
pdf(paste(type,"is_fig1.pdf",sep="_"))
matplot(t(t(insertCounts[1:250,])/colSums(insertCounts)),type="l",lty=1,lwd=3,col=colors,ylab="",main=type)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()

