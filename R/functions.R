library(wesanderson)

check_infile=function(file){
  ## sort so that input last, rep 1 before 2 and samples in order 
  file_sample=split(file,file$V3)$sample
  file_input=split(file,file$V3)$input
  file_sample$V3
  if (sum(order(file_sample$V4)==1:nrow(file_sample))!=nrow(file_sample)){
      print("types not ordered in infile, reordering")
      file_sample = file_sample[order(file_sample$V4),]
  }
  if (sum(order(file_sample$V4,file_sample$V5)==1:nrow(file_sample))!=nrow(file_sample)){
      print("replicates not ordered in infile, reordering")
      file_sample = file_sample[order(file_sample$V4,file_sample$V5),]
  }
  new_file = rbind(file_sample,file_input) 
  new_file
}




defineColors=function(file){
  file = check_infile(file)
  
  two_times_two_group_colors = c(wes_palette(n=4, name="Zissou")[c(2:1,3:4)]) ## could add third 
  two_input = c("darkblue","orange")
  one_input = c("darkgrey")
  genericColors = list(c("lightblue","darkblue"),c("lightgreen","darkgreen"),c("pink","darkred"),c("yellow","darkorange"))
print(head(file))  
  sp =  lapply(lapply(split(file,file$V3),droplevels),function(x) split(x,x$V4))
	print(sp)  
file_sample=droplevels(split(file,file$V3)$sample)
  #file_input=droplevels(split(file,file$V3)$input)
  nrCond = length(sp$sample)
  #nrInp = length(sp$inputi)
nrInp=0
  nrReps = lapply(sp,function(x) lapply(x,nrow))
  if(nrCond<3 & max(unlist(nrReps))<3 & nrInp<3){ ## use wes colors (less than 3 conditions, less than 3 replicates, max 2 input)
    file_sample$col= two_times_two_group_colors[c(1:nrReps$sample[[1]],c(3:(2+nrReps$sample[[2]])))]
    #if(nrInp==2 & max(unlist(nrReps$input))==1){
     # file_input$col=two_input
   #   } else
      # file_input$col=one_input
  } else {
      cc=list()
      for (i in 1:length(nrReps$sample))
        cc[[i]]=colorRampPalette(genericColors[[i]])(nrReps$sample[[i]])
      file_sample$col = unlist(cc)
       #file_input$col=one_input
  }
  #new = rbind(file_sample,file_input)
new = file_sample  
new
  }
    
#input = c("SCN_rep2.marked_duplicates.bam.numbers.tab", "SCN_rep1.marked_duplicates.bam.numbers.tab" ,"VCN_rep1.marked_duplicates.bam.numbers.tab" , "VCN_rep2.marked_duplicates.bam.numbers.tab",     
          # "VCN_input.marked_duplicates.bam.numbers.tab")
#file =read.table("files.tab")
checkInput=function(input,file){
  file = check_infile(file)
  fileName = file[,1]
  if (sum(sapply(fileName,function(x) grep(x,input))==1:length(input))!=length(input)){
    print("reordering")
    ord =  sapply(fileName,function(x) grep(x,input))
    input = input[ord]
  }
  input
}

findRepCondOrder=function(files,fileF){
  samF = droplevels(subset(fileF,fileF$V3=="sample"))
  gr=as.numeric(samF$V4)
  rep=samF$V5
  return(data.frame(group=gr,rep=rep))
}

## replicate venn diagrams across 2 samples
plotReplicateVenn=function(sumList,v_colors,gr,group,cat){
  if (sum(gr==group)<2 | sum(gr==group)>3){
    stop("has to be exactly 2 or 3 samples in the pairwise")
  }
  tab = table(sumList[[2]][[group]])
  if(sum(gr==group)==2)
    vp = draw.pairwise.venn(tab["1.0"]+tab["1.1"],tab["0.1"]+tab["1.1"],tab["1.1"],fill=v_colors[which(gr==group)] ,category = cat[which(gr==group)])
  if(sum(gr==group)==3)
    vp = draw.triple.venn(area1=tab["1.0.0"]+tab["1.1.0"]+tab["1.1.1"]+tab["1.0.1"],
                            area2=tab["0.1.0"]+tab["1.1.0"]+ tab["0.1.1"]+tab["1.1.1"],
                            area3=tab["0.0.1"]+tab["1.0.1"]+ tab["0.1.1"]+tab["1.1.1"],
                            n12=tab["1.1.0"]+tab["1.1.1"],
                            n23=tab["0.1.1"]+tab["1.1.1"],
                            n13=tab["1.0.1"]+tab["1.1.1"],
                            n123 = tab["1.1.1"],
                            fill=v_colors[which(gr==group)] , category = cat[which(gr==group)])
  grid.draw(vp)
}  

## approach 1 common if in all, unique if in all rep (1 or zero other)

defineA = function(overlaps,gr,m){
  p=list()
  for (i in 1:max(gr)){
    p[[i]]=rowSums(overlaps[,gr==i,drop=FALSE])>=m[i]  
  }
  p
}

plotAppVenn=function(p,o_col,cat){
  if(length(p)>3 | length(p)<2){
    stop("min2, max 3")
  }
  ps = do.call("cbind",p)
  if(length(p)==2){
    n12 = sum(rowSums(ps)==2)
    a1 = sum(ps[,1]==TRUE)
    a2 = sum(ps[,2]==TRUE)
    vp = draw.pairwise.venn(a1,a2,n12,fill=o_col[1:2],category=cat)
    grid.draw(vp)
    }
  if(length(p)==3){
    n123 =  sum(rowSums(ps)==3)
    a1 = sum(ps[,1]==TRUE)
    a2 = sum(ps[,2]==TRUE)
    a3 = sum(ps[,3]==TRUE)
    n12 = sum(ps[,1]+ps[,2]==2)
    n23 = sum(ps[,3]+ps[,2]==2)
    n13 = sum(ps[,1]+ps[,3]==2)
    vp = draw.triple.venn(a1,a2,a3,n12,n23,n13,n123,fill=o_col,category=cat)
    grid.draw(vp)
    }
}

