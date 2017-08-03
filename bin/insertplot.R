args = commandArgs(trailingOnly=TRUE)


## first arg is type

type=args[1]

## rest files

files = args[2:length(args)]

################
## distribution fit function (should be fitted on fragment sizes (or my eqvivalent))

my_gamma_fit=function(X,o,p1,p2,p3){
  k=p1
  theta=p2
  a=p3
  x_mod=X-o
  #print(x_mod)
  res=matrix(0,nrow=length(x_mod))
  if (k>1){
    nz = x_mod[x_mod>=0]
    nz = which(x_mod>=0)
  }
  else{
    nz = x_mod[x_mod>0]
    nz = which(x_mod>0)
  }
  res[nz]=a*x_mod[nz]^(k-1)*exp(-x_mod[nz]/theta) / (theta^k*gamma(k))
  res
}


run_python_opt = function(y){
  unlink('y_data.txt')
  write.table(y,file="y_data.txt",col.names = FALSE, row.names = FALSE)
  opt = system2("python", args="../../../bin/rPython.py",stdout = TRUE)
  if(length(opt)>1){
    print(opt[1])
    opt = opt[length(opt)]
  }
  opt = as.numeric(gsub("\\]","",gsub("\\[","",strsplit(opt,",")[[1]])))
  names(opt)=c("p1","p2","p3","o")
  opt
}



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
print(dim(insertCounts))

print(dim(is[[1]]))
print(dim(is[[2]]))

for (i in names(is)){
print(i)
print(rownames(is[[i]]))
print(unlist(is[[i]])
  insertCounts[rownames(is[[i]]),i]=unlist(is[[i]])
}

print(t(t(insertCounts[1:250,])/colSums(insertCounts)))

#create plots
pdf(paste(type,"is_fig1.pdf",sep="_"))
matplot(t(t(insertCounts[1:250,])/colSums(insertCounts)),type="l",lty=1,lwd=3,col=colors,ylab="",main=type)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()


CairoPNG(paste(type,"is_fig1.png",sep="_"))
matplot(t(t(insertCounts[1:250,])/colSums(insertCounts)),type="l",lty=1,lwd=3,col=colors,ylab="",main=type)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()


para_test = list()
test_fit0 = matrix(0,ncol=length(files),nrow=250)
colnames(test_fit0)=names(is)
for (i in files){
  para_test[[i]] = run_python_opt(insertCounts[36:115,i]/sum(insertCounts[,i])) ## 36:115 is what is used in NucleATAC
  test_fit0[,i] = my_gamma_fit(X=matrix(1:250),o=para_test[[i]]["o"],para_test[[i]]["p1"],para_test[[i]]["p2"],para_test[[i]]["p3"])
  }

pdf(paste(type,"is_fig2.pdf",sep="_"))
matplot(test_fit0,type="l",lty=1,col=colors,lwd=2,main=type)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()

CairoPNG(paste(type,"is_fig2.png",sep="_"))
matplot(test_fit0,type="l",lty=1,col=colors,lwd=2,main=type)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()



