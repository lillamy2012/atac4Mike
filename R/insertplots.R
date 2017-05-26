args = commandArgs(trailingOnly=TRUE)

source("../R/functions.R")
if (length(args)==0){
  stop("Provide tab filewith files to read")
} else if (length(args)>1) {
  print ("more than one input file, only first will be processed")
} else
  input=args[1]



### load libraries
library(wesanderson)
library(Cairo)

### define colors
#colors = wes_palette(n=4, name="Zissou")[c(2:1,3:4)]
#colors = c(wes_palette(n=4, name="Zissou")[c(2:1,3:4)],"darkblue","orange") ## need to change
### in argument with files to use
#input = "../tables/R_uniq_filtered.bam.IS.tab"
files=read.table(input)
fileF=read.table("../files.tab",comment.char="")
colors=defineColors(fileF)
files  = files[,1]
files = checkInput(input=files,file=colors)
print(colors)
print(files)
colors = colors$col

###
name = basename(input)
name = sub("R_","",name)
lname = strsplit(name,"\\.")[[1]]
exclude = c("bam","tab","txt")
new = paste(lname[!lname%in%exclude],collapse = "_")


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
  opt = system2("python", args="../python/rPython.py",stdout = TRUE)
  if(length(opt)>1){
    print(opt[1])
    opt = opt[length(opt)]
  }
  opt = as.numeric(gsub("\\]","",gsub("\\[","",strsplit(opt,",")[[1]])))
  names(opt)=c("p1","p2","p3","o")
  opt
}

### read in data

is=list()
for (i in files){
  is[[i]] = read.table(i)
  is[[i]]$V1 = is[[i]]$V1-8 ## 
  rownames(is[[i]])=is[[i]]$V1
  is[[i]]=is[[i]][,2,drop=FALSE]
}

## to matrix
nrows = max(sapply(is,function(x)  max(as.numeric(rownames(x)))))
ncols = length(is)
insertCounts = matrix(0,ncol=ncols,nrow=nrows)
colnames(insertCounts)=names(is)
rownames(insertCounts)=1:nrows

for (i in names(is)){
  insertCounts[rownames(is[[i]]),i]=unlist(is[[i]])
}
 
#create plots
pdf(paste(new,"fig1.pdf",sep="_"))
matplot(t(t(insertCounts[1:250,])/colSums(insertCounts)),type="l",lty=1,lwd=3,col=colors,ylab="",main=new)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()

CairoPNG(paste(new,"fig1.png",sep="_"))
matplot(t(t(insertCounts[1:250,])/colSums(insertCounts)),type="l",lty=1,lwd=3,col=colors,ylab="",main=new)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()


para_test = list()
test_fit0 = matrix(0,ncol=length(files),nrow=250)
colnames(test_fit0)=names(is)
for (i in files){
  para_test[[i]] = run_python_opt(insertCounts[36:115,i]/sum(insertCounts[,i])) ## 36:115 is what is used in NucleATAC
  test_fit0[,i] = my_gamma_fit(X=matrix(1:250),o=para_test[[i]]["o"],para_test[[i]]["p1"],para_test[[i]]["p2"],para_test[[i]]["p3"])
  }

pdf(paste(new,"fig2.pdf",sep="_"))
matplot(test_fit0,type="l",lty=1,col=colors,lwd=2,main=new)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()

CairoPNG(paste(new,"fig2.png",sep="_"))
matplot(test_fit0,type="l",lty=1,col=colors,lwd=2,main=new)
legend("topright",legend = files, lwd=3,col=colors)
dev.off()

