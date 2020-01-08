# GENCODE id conversion file. Contains the 'attr' variable.
load("/data/genomes/gene_annotations_v19.Rdata")

# Create filter of X chromosome genes. 
f.x = attr$chr == "chrX"


# Load SNP counts data (eg.). See other scripts that generated these. Note, new GATK output means some columns are modified. Will specify. 
load("EGAD00001001086.distributions.chrX.Rdata")


# Calculate skew ratios (ref and alt)
list.skew = list()
for( j in 1:length(list.vcf)) {
  if( is.null(list.vcf[[j]]) ) { next }
  f1 = rowSums(list.vcf[[j]][,2:3] > 0 ) == 2

  skews = list.vcf[[j]][f1,1:6]
  skews[,4] = rowSums(list.vcf[[j]][f1,2:3])
  skews[,5:6] = list.vcf[[j]][f1,2:3]

  list.skew[[j]] = skews
}


Ns = sapply(1:length(list.skew), function(i) dim(list.skew[[i]])[1] )
A = sapply(1:length(list.skew), function(i) is.null(Ns[[i]] ) )  * 1
for( i in 1:length(Ns) ) { if( is.null(Ns[[i]]))  { A[i] = 0 } else { A[i] = Ns[[i]] }}
Ns = A

# Find SNP/gene overlaps. Note, takes the first gene if there are genes on positive and negative strands. Need to fix/select. 
list.skew2 = list.skew # save the original skew - can delete this
for( j in 1:length(list.skew)) {
  if( is.null(list.skew[[j]]) ) { next }

  genes.snps = sapply( 1:dim(list.skew[[j]])[1], function(i) which(list.skew[[j]][i,1] > attr[f.x,2]  &  list.skew[[j]][i,1] < attr[f.x,3] ))
  for(i in 1:length(genes.snps)) { if( length(genes.snps[[i]])==0 ) { genes.snps[[i]] = NA } ;  if( length(genes.snps[[i]]) > 1  ) { genes.snps[[i]] = genes.snps[[i]][1]  }   }
  list.skew[[j]] = cbind(list.skew[[j]], unlist(genes.snps), attr[f.x,][unlist(genes.snps),]  )
}


# save(  list.skew, Ns, file="skew.est.filt.Rdata")





## Pick max powered SNPs
list.skew.max = list()
for(j in 1:length(list.skew)){ 
  if(Ns[j] > 0 ){  
    a = list.skew[[j]][,4]
    b = list.skew[[j]][,13]
    
    abi = tapply(a,b, which.max)
    abi = abi[!is.na(abi)]
    abi = cbind(names(abi), abi )
    maxtest  = lapply(1:dim(abi)[1], function(i) list.skew[[j]][ which(list.skew[[j]][,13]==abi[i,1]),][ as.numeric(abi[i,2] ),] )

    list.skew.max[[j]] = do.call(rbind, maxtest)
  }
}
 rm(maxtest)

## Store as a separate matrices/lists
ratios.max.genes =  list() 
ratios.max = list()
for(j in 1:length(list.skew.max)){ 
  if(Ns[j] > 0 ){
    ratios.max.genes[[j]] = list.skew.max[[j]][,c(13:17,8:12)]
    ratios.max[[j]] = list.skew.max[[j]][,5:6]
  }
}

 save(list.skew.max, ratios.max.genes,ratios.max, Ns, file="skew.est.max.genes.Rdata" )


##############################################
### Or you can load a few of the files from here. 
##############################################
# eg. 
load("GTEx.blood.skew.est.max.genes.Rdata")


# Functions to fold and unfold ratios 
folded <- function(x) { apply( cbind(x,1-x), 1, max)} 
unfold <- function(x) { c(x,1-x) } 

# Maximum likelihood estimate function. Note, x is a global variable! Should fix but this lets me use it in the bootstrap later.  
mle_folded <- function(){ 
  mus = seq(0.5,1, by = 0.001)
  sigmas = seq(0,0.5, by = 0.01)
  
  mles = sapply(mus, function(mu)
    sapply(sigmas, function(sigma)
      -sum( log(dnorm( x,   mean = mu, sd = sigma )
                + dnorm( x,   mean = 1-mu, sd = sigma ) )  )))
  mles[!is.finite(mles)] = NA
  coefs = which(mles==min(mles, na.rm=T), arr.ind=T)
  return (list( mus[coefs[2]] , sigmas[coefs[1]]))
} 


# Estimate skews 
# Probabilities to test   
p = 50:100/100

fittfold.list = list() 
est_fold.list = list() 

# Note, mle_folded() needs x which is the list of skew ratios. 
for(j in 1:length(list.skew.max)){ 
  fittfold.list[[j]] = list()
  est_fold.list[[j]] = list() 

  if(Ns[j] > 0) { 
    for(i in 1:2){ 
      x = folded(ratios.max[[j]][,i])
      
      fittfold.list[[j]][[i]] =  mle_folded() 
      est_fold.list[[j]][[i]] = p[which.max( dnorm(p, mean = fittfold.list[[j]][[i]][[1]] , sd = fittfold.list[[j]][[i]][[2]] )+ 
                                               dnorm(p, mean = 1-fittfold.list[[j]][[i]][[1]] , sd = fittfold.list[[j]][[i]][[2]] ) )]
      
    }   
  }
} 

library(viridis)

# Note, this does both ref and alt which should be the same. Was just a sanity check.  
est_fold  = sapply(which(Ns>0), function(j) sapply(1:2, function(i) est_fold.list[[j]][[i]][[1]])) 

# Plot distribution (folded)
hist(est_fold[1,], xlab="Skew",  main="Folded normal estimates", col=magma(10)[2], border=NA)

# Plot distribution (unfolded)
hist( unfold(est_fold[1,]), xlab="Skew",  main="Folded normal estimates", col=magma(10)[5], border=NA, freq=F, xlim=c(0,1))

# Plot curve estimates 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/2) ) ,   lwd=2, col=magma(10)[1] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/4) ) ,   lwd=2, col=magma(10)[2] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/8)) ,    lwd=2, col=magma(10)[3] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/16) ) ,  lwd=2, col=magma(10)[4] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/32) ) ,  lwd=2, col=magma(10)[5] ) 


# EStimate variance of distribution 
var = sd((unfold(est_fold[1, ])) , na.rm=T ) ^ 2
# Estimate cell count based on Ncells = pq/var. Since p == q == 0.5, Ncells = 0.5 * 0.5 / var, Ncells = 0.25/var  
Ncells = 0.25/var

# Bootstrap 
library(boot)
foo <- function(i,j) {
  return( 0.25/(sd(unfold(i[j]))^ 2))
} 
myBootstrap <- boot(array(est_fold[1,]), foo, R=1000 )
bootCI = boot.ci(myBootstrap, index=1)


pi = 0.5
var_est = (pi * (1-pi)) / (1:250)

# Plot curve estimate
plot( (var_est), pch=19, xlab="N cells", ylab="Var")
abline(h=var, col=2, lwd=2)
abline(v=Ncells, col=2,lwd=2)
polygon(c(bootCI$percent[4],bootCI$percent[5],bootCI$percent[5],bootCI$percent[4]), c(0,0,1,1) , 
        col=makeTransparent(magma(5)[4]),
        border=NA)



