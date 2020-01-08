library(plyr)

count_pairs <- function( data){
  freqs = count(data)
  pairs = freqs[,1] > 0 &  freqs[,2] > 0 &  freqs[,3] > 1
  o = order(freqs[pairs,3], decreasing = T)
  return(freqs[pairs,][o,]) 
}

# Carrier 1 data
#load("/data/XCGD/SCGC-GILL-JG-03_2/carrier1.temp_by_cell.Rdata")
# Carrier 2 data
#load("/data/XCGD/SCGC-GILL-JG-04_2/carrier2.temp_by_cell.Rdata")

update_freqs <- function(mat ){
  
  A = mat
  B = mat 
  
  A[A>0] = 1
  num2 = t(A) %*% (A)
  diag(num2) = 0
  num2[num2==0] = NA
  a = (num2[row(num2) > col(num2)] )
  a = a[!is.na(a)]
  freq = count(a)
  o = order(freq[,1], decreasing = T)
  freq = freq[o,]
  hist(a[a>1], breaks=100, xlab="Number of cells with common pair of SNPs", main="Carrier 1")
  #  hist(a[a>19], breaks=100, xlab="Number of cells with common pair of SNPs", main="Carrier 1")
  nn = which(freq[,1]==10)
  ind.all = do.call(rbind, lapply( freq[1:nn,1], function(top) which(num2==top, arr.ind=T)) )
  ind.all = unique(t(apply(ind.all,1,sort)) )
  
  ind = ind.all[1,]
  return(ind)
}


update_freqs2 <- function(mat ){
  
  A = mat
  B = mat 
  
  A[A>0] = 1
  num2 = t(A) %*% (A)
  diag(num2) = 0
  num2[num2==0] = NA
  a = (num2[row(num2) > col(num2)] )
  a = a[!is.na(a)]
  freq = count(a)
  o = order(freq[,1], decreasing = T)
  freq = freq[o,]
  hist(a[a>1], breaks=100, xlab="Number of cells with common pair of SNPs", main="Carrier 1")
  # hist(a[a>19], breaks=100, xlab="Number of cells with common pair of SNPs", main="Carrier 1")
  nn = which(freq[,1]==10)
  ind.all = do.call(rbind, lapply( freq[1:nn,1], function(top) which(num2==top, arr.ind=T)) )
  ind.all = unique(t(apply(ind.all,1,sort)) )
  
  #ind = ind.all[1,]
  return(ind.all)
}
count_pairs <- function( data){
  freqs = count(data)
  pairs = freqs[,1] > 0 &  freqs[,2] > 0 &  freqs[,3] > 1    
  o = order(freqs[pairs,3], decreasing = T)
  return(freqs[pairs,][o,]) 
}

count_pairs_pvals <- function( data){
  freqs = count(data)
  pairs = freqs[,1] > 0 &  freqs[,2] > 0  
  freqs = freqs[pairs,]
  
  As = unique(freqs[,1])
  Bs = unique(freqs[,2])
  table = matrix(0, ncol=length(As), nrow=length(Bs))
  colnames(table) = As
  rownames(table) = Bs
  if( length(As) == 2 & length(Bs) ==  2 ){
  for( i in 1:dim(freqs)[1]){
    Ai = freqs[i,1]
    Bi = freqs[i,2]
    Ci = freqs[i,3]
    
    ii = which(As == Ai)
    ij = which(Bs == Bi)
    table[ij,ii] = Ci 
  }
  
  pval = fisher.test(table)$p.val
  return(list(table,pval)) 
 }
 return(list(NA,NA))
}

A = cell_by_snp.hap.num[,filt]
filtr = (rowSums(A > 0 ) > 1  )
A = cell_by_snp.hap.num[filtr,filt]
B = cell_by_snp.hap.num[filtr,filt]
A = A[,filt2]
B = B[,filt2]


A[A>0] = 1
num = t(A) %*% (A)
num2 = num
diag(num2) = 0
num2[num2==0] = NA
a = (num2[row(num2) > col(num2)] )
a = a[!is.na(a)]
freq = count(a)
o = order(freq[,1], decreasing = T)
freq = freq[o,]
hist(a[a>1], breaks=100, xlab="Number of cells with common pair of SNPs", main="Carrier 1")
hist(a[a>19], breaks=100, xlab="Number of cells with common pair of SNPs", main="Carrier 1")
nn = which(freq[,1]==10)
ind.all = do.call(rbind, lapply( freq[1:nn,1], function(top) which(num2==top, arr.ind=T)) )
ind.all = unique(t(apply(ind.all,1,sort)) )

dim(ind.all)
snpss = list()
mat.build.merge = B 
pj = list() 
hap_i = c(6,6+1e6) 
ind = ind.all[1,]
k = 1 

fishp = sapply(1:dim(ind.all)[1], function(i) count_pairs_pvals(mat.build.merge[,ind.all[i,]])[[2]] )
np = sapply(1:dim(ind.all)[1], function(i) sum(count_pairs_pvals(mat.build.merge[,ind.all[i,]])[[1]]) )
#while( length(ind) > 0  ){ 
o = order(fishp)

  #ind = ind.all[o[1],]
  ind = ind.all[o[2],]
  
  pj[[k]] = count_pairs(mat.build.merge[,ind]) 
  pj[[k]]
  ji = ind[1]
  jj = ind[2]
  n = dim(pj[[k]])[1]
  
  if( n>0 ){ 
    for(i in 1:2){
      fji = mat.build.merge[,ji] == as.numeric(pj[[k]][i,1]) 
      fjj = mat.build.merge[,jj] == as.numeric(pj[[k]][i,2])
      
      fji.miss = which(mat.build.merge[fji,jj] == 0) 
      if( length(fji.miss) >0){
        mat.build.merge[fji,jj][fji.miss] = hap_i[i]
      }
      fjj.miss = which(mat.build.merge[fjj,ji] == 0)
      if( length(fjj.miss) > 0){
        mat.build.merge[fjj,ji][fjj.miss] = hap_i[i]
      }
      mat.build.merge[fjj,jj] = hap_i[i]
      
      hap_i[i] = hap_i[i] + 1 
      
    } 
    colnames(mat.build.merge)[jj] = paste(colnames(mat.build.merge)[c(ji,jj)], collapse="_")
    mat.build.merge = mat.build.merge[,-ji]
  }
  k = k + 1 
  ind.all = update_freqs2(mat.build.merge)
  fishp = sapply(1:dim(ind.all)[1], function(i) count_pairs_pvals(mat.build.merge[,ind.all[i,]])[[2]] )
  np = sapply(1:dim(ind.all)[1], function(i) sum(count_pairs_pvals(mat.build.merge[,ind.all[i,]])[[1]]) )
  
  plot(np, -log10(fishp), pch=19, xlab="Number of cells (minimum 10)", ylab="-log10 pvalue")
  # abline(0,1, col=2)
  #count_pairs(mat.build.merge[,ind]) 
  
#}


