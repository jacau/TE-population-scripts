#this is an updated version of the old file which is from 2015

dist_to_gene <- function(gchr,tchr){
  
  for(t in 1:nrow(tchr)){
    
    te_st = tchr$start[t]
    te_e = tchr$end[t]
    len = te_e - te_st
    
    if(te_st < gchr$start[1]){
      
      if(len < 75){
        tchr$gdist[t] = max(gchr$start[1] - te_st,0)
      }else{
        tchr$gdist[t] = max(gchr$start[1] - te_e,0)
      }
      next
      
    }else if(te_st > gchr$start[nrow(gchr)]){
      tchr$gdist[t] = max(te_st - gchr$end[nrow(gchr)],0)
      next
    }
    
    temp = rbind(gchr,c(tchr$scaff[t],'TE',te_st,te_e))
    temp$start = as.numeric(temp$start)
    temp$end = as.numeric(temp$end)
    temp = temp[order(temp$start),]
    indx = which(temp$type == 'TE')
    
    gene_end = temp$end[indx-1]
    
    if(indx-2 > 0){
      #in some cases (because of birectionality of genes) after sort by start, end of previous gene overhangs end of next gene
      if(temp$end[indx-2] > temp$end[indx-1]){ 
        gene_end = temp$end[indx-2]
      }
    }
    
    if(len < 75){ #is a non-reference TE
      
      d1 = te_st - gene_end
      d2 = temp$start[indx+1] - te_st
      
      if(d1 <= 0 | d2 <= 0){ #TE is in the gene
        dist = 0
      }else{
        dist = min(d1,d2)
      }
      
    }else{ #is a reference TE (only difference is to use t_e in d2 instead of te_st)
      
      d1 = te_st - gene_end
      d2 = temp$start[indx+1] - te_e
      
      if(d1 <= 0 | d2 <= 0){ #TE is in the gene
        dist = 0
      }else{
        dist = min(d1,d2)
      }
      
    }
    
    tchr$gdist[t] = dist
    
  }
  return(tchr)
  
}

#command line arguments: te_db file, gene location file, output file name
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 3){
  stop("Rscript dist_to_gene.R [TE file] [gene annotation] [output file]", call. = FALSE)
}


te = read.table(args[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
gene = read.table(args[2], header=TRUE, sep="\t", stringsAsFactors=FALSE)

gene$Direction = NULL
gene$PacID = NULL
names(gene) = c('scaff','type','start','end')

te$gdist = 0

scaff_gene = c('scaffold_1','scaffold_2','scaffold_3','scaffold_4','scaffold_5','scaffold_6','scaffold_7','scaffold_8')
scaff_te = c('S1','S2','S3','S4','S5','S6','S7','S8')

for(x in 1:8){
  tchr = subset(te,te$scaff == scaff_te[x])
  gchr = subset(gene,gene$scaff == scaff_gene[x])
  gchr = gchr[order(gchr$start),]
  temp = dist_to_gene(gchr,tchr)
  if(!x==1){
    df = rbind(df,temp)
  }else{
    df = temp
  }
}

write.table(df, args[3], sep='\t', row.names = FALSE)
