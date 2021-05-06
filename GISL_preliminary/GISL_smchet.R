
############################# FUNCTIONS ################################

subclone.dirichlet.gibbs <- function(C = 20, a, d, copyNumberCoef = rep(2, length(a)), iter = 1000) {
  # C: maximum number of subclones 
  # a: a vector of variant read counts
  # d: a vector of total read counts
  # iter: number of iteration of Gibbs sampler
  
  N <- length(a)
  print(paste("num.muts=", N, sep=""))
  
  # Hyperparameters for alpha
  A <- B <- 0.01
  
  # set up data formats for recording iterations
  alpha <- rep(NA, iter)
  pi.h <- matrix(NA, nrow = iter, ncol = C)     # {Phi_k} for k = 1,...,K, i.e., cellular proportion of each cluster
  V.h <- matrix(NA, nrow = iter, ncol = C)      # stick-breaking
  z.i <- matrix(NA, nrow = iter, ncol = N)      # data assignments
  Pr.z <- matrix(NA, nrow = N, ncol = C)        # Pr(z.n = k), used for sample {zn}
  mutBurden <- array(NA, dim = c(iter, C, N))   # expected fraction of reads containing mutation, used as p in Binmomial(n,p)
  
  # randomize starting positions
  alpha[1] <- 1
  V.h[1,] <- c(rep(0.5, C-1), 1)  # Truncated Dirichlet Process(TDP): last element is set to 1 so that pi.h[,k] = 0 for all k>C
  V.h[1:iter, C] <- rep(1, iter)
  z.i[1,] <- sample.int(C, size = N, replace = TRUE)
  
  lower <- min(copyNumberCoef*a/d)
  upper <- max(copyNumberCoef*a/d)
  difference <- upper - lower
  pi.h[1,] <- runif(C, lower<-(lower-difference/10), upper<-(upper+difference/10))
  for(c in 1:C){
    tmp <- pi.h[1,c]/copyNumberCoef
    tmp[which(tmp<0.000001)] <- 0.000001
    tmp[which(tmp>0.999999)] <- 0.999999
    mutBurden[1,c,] <- tmp
    # mutBurden[1,c,] <- pi.h[1,c]/copyNumberCoef
    # mutBurden[1,c, mutBurden[1,c,]>0.499999] <- 0.499999
  }
  
  # update cluster assignment for each mutation
  for(i in 2:iter) {
    if(i %% 100 == 0) print(i)
    
    # independently sample zn from Pr(zn=k)
    for(n in 1:N) {
      # use log-space
      Pr.z[n,1] <- log(V.h[i-1,1]) + a[n]*log(mutBurden[i-1,1,n]) + (d[n]-a[n])*log(1-mutBurden[i-1,1,n])
      # Pr.z[n,1] <- log(V.h[i-1,1]) + dbinom(a[n],d[n],prob = mutBurden[i-1,1,n], log = TRUE)
      Pr.z[n,2:C] <- sapply(2:C, function(j) {log(V.h[i-1,j])+sum(log(1-V.h[i-1, 1:(j-1)])) + a[n]*log(mutBurden[i-1,j,n]) + (d[n]-a[n])*log(1-mutBurden[i-1,j,n])})
      # Pr.z[n,2:C] <- sapply(2:C, function(j) {log(V.h[i-1,j])+sum(log(1-V.h[i-1, 1:(j-1)])) + dbinom(a[n],d[n],prob = mutBurden[i-1,j,n], log = TRUE)})
      # library("matrixStats)
      Pr.z[n,] <- exp(Pr.z[n,] - logSumExp(Pr.z[n,],na.rm = TRUE))  
      Pr.z[n,is.na(Pr.z[n,])] = 0
    }
    z.i[i,] <- sapply(1:N, function(j) {sum(rmultinom(1,1,Pr.z[j,]) * (1:C))})
    
    # independently sample V.h
    V.h[i,1:(C-1)] <- sapply(1:(C-1), function(j) {rbeta(1, 1+sum(z.i[i,] == j), alpha[i-1]+sum(z.i[i,] > j))})
    V.h[i,which(V.h[i,1:(C-1)] == 1)] <- 0.999    # prevent one stick from taking all the weight
    
    # independently sample pi.h from the posterior distribution in the same family as the base distribution
    pi.h[i,] <- runif(C, lower, upper)
    mutBurden[i,,] <- mutBurden[i-1,,]
    for(c in unique(z.i[i,])) {
      # the mean of gamma distribution is shape/rate.
      pi.h[i,c] <- rgamma(1, shape = sum(a[z.i[i,] == c]), rate = sum((d/copyNumberCoef)[z.i[i,] == c]))
      tmp <- pi.h[i,c]/copyNumberCoef
      tmp[which(tmp<0.000001)] <- 0.000001
      tmp[which(tmp>0.999999)] <- 0.999999
      mutBurden[i,c,] <- tmp
      # mutBurden[i,c,] <- pi.h[i,c]/copyNumberCoef
      # mutBurden[i,c, mutBurden[i,c,]>0.499999] <- 0.499999
    }
    
    # update alpha
    alpha[i] <- rgamma(1, shape = A+C-1, rate = B-sum(log(1-V.h[i,1:(C-1)])))
    
  }
  return(list(z.i=z.i, V.h = V.h, Phi.h = pi.h, alpha = alpha, a = a, d = d))
}

Gibbs.subclone.density.est <- function(GS.data, density.file = "density.txt", density.smooth = 0.1,
                                       burn.in = 300, y.max = 5){
  # GS.data is the list output from function -- subclone.dirichlet.gibbs()
  # density.smooth is the parameter used in R's density() function
  
  V.h <- GS.data$V.h
  Phi.h <- GS.data$Phi.h
  no.iter <- nrow(V.h)
  wts <- matrix(NA, nrow = no.iter, ncol = ncol(V.h))
  wts[,1] <- V.h[,1]
  wts[,2] <- V.h[,2] * (1-V.h[,1])
  for (i in 3:dim(wts)[2]) 
    wts[,i] <- apply(1-V.h[,1:(i-1)], MARGIN=1, FUN=prod) * V.h[,i]
  
  xx <- density(c(Phi.h[burn.in-1,]), weights=c(wts[burn.in,]) / sum(c(wts[burn.in,])), adjust=density.smooth, from=0, to=1)$x
  
  post.density <- matrix(NA, ncol = no.iter-burn.in+1, nrow = 512) # nrow is equal to the output of function density()
  for(i in burn.in:no.iter)
    post.density[,i - burn.in + 1] <- density(c(Phi.h[i-1,]), weights=c(wts[i,]) / sum(c(wts[i,])), adjust=density.smooth, from=0, to=1)$y
  
  # take the median of density
  yy <- apply(post.density, MARGIN = 1, FUN = quantile, probs = 0.5)
  write.table(cbind(xx,yy),density.file,sep="\t",col.names=c("mutation.burden","median.density"),row.names=F,quote=F)
}


getClusterAssignment <- function(GS.data, density.file, window.size = 20, burn.in = 300) {
  mutReads <- GS.data$a
  totalReads <- GS.data$d
  no.muts = length(mutReads)
  
  z.i = GS.data$z.i
  V.h = GS.data$V.h
  Phi.h = GS.data$Phi.h
  iter = nrow(z.i)
  
  density = read.table(density.file, row.names = NULL, header = TRUE, sep = "\t")
  
  # use sliding window to find local optima
  localOptima = NULL
  peak.indices = NULL
  for(i in (1+window.size):(nrow(density)-window.size)) {
    if(density$median.density[i] == max(density$median.density[(i-window.size):(i+window.size)])) {
      localOptima = c(localOptima, density[i,1])
      peak.indices = c(peak.indices, i)
    }
  }
  
  # assign mutations to clusters
  no.optima = length(localOptima)
  if(no.optima>1){
    boundary = array(NA, no.optima-1)
    for(i in 1:(no.optima-1)) {
      min.density = min(density$median.density[(peak.indices[i]+1):(peak.indices[i+1]-1)])
      min.indices = intersect(which(density$median.density == min.density), (peak.indices[i]+1):(peak.indices[i+1]-1))
      # positions of each minimum between pairs of optima
      boundary[i] = (density[max(min.indices),1] + density[min(min.indices),1])/2
    }
    
    sampledIters = (burn.in+1):iter
    # do not use the initial state
    sampledIters = sampledIters[sampledIters != 1]
    if(length(sampledIters) > 1000) 
      sampledIters = floor(burn.in + (1:1000)*(iter - burn.in) / 1000)
    
    z.i = data.matrix(z.i)
    mutation.preferences = array(0, c(no.muts, no.optima))
    for(s in sampledIters) {
      temp.preferences = array(0, c(no.muts, no.optima))
      for(c in unique(z.i[s,])){
        # use the relative position between Phi.h and boundary to determine the cluster it belongs to
        bestOptimum = sum(Phi.h[s,c] > boundary)+1
        temp.preferences[z.i[s,]==c, bestOptimum] = temp.preferences[z.i[s,]==c, bestOptimum] + 1
      }
      iter.preferences = t(sapply(1:no.muts, function(p,i) {as.integer(p[i,]==max(p[i,]))/sum(p[i,]==max(p[i,]))}, p = temp.preferences))
      mutation.preferences = mutation.preferences + iter.preferences
    }
    mutation.preferences = mutation.preferences/length(sampledIters)
    most.likely.cluster = sapply(1:no.muts, function(i) which.max(mutation.preferences[i,]))
  } else {
    most.likely.cluster = rep(1, no.muts)
  }
  
  return(list(localOptima = localOptima, most.likely.cluster = most.likely.cluster))
  
}


### MAIN CODE ####

args<-commandArgs(TRUE)
suppressMessages(require("matrixStats"))  # suppress message to avoid error when running on Galaxy

#################################### Subchallenge 1 & 2 ################################################

ssmdat = read.table(args[1],sep='\t',comment.char='#', stringsAsFactors = FALSE)
colnames(ssmdat) = c("CHROM","POS","ID","REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT","normal", "tumor")
tumor_stat = do.call(rbind, strsplit(as.vector(ssmdat[,"tumor"]), split = ":", fixed = TRUE))
colnames(tumor_stat) = strsplit(as.vector(unique(ssmdat[, "FORMAT"])), ":")[[1]]
rownames(ssmdat) = paste(ssmdat[,"CHROM"], ssmdat[,"POS"], sep = "_")

cnvdat <- read.table(args[2], header = TRUE, stringsAsFactors = FALSE)
rownames(cnvdat) = paste(cnvdat[,1], cnvdat[,2], cnvdat[,3], sep = "_")

# use battenberg data to estimate cellularity0. NOTE: not the final output
id.xy <- c(which(cnvdat$chr == "X"),which(cnvdat$chr == "Y"))
id.sub <- which(cnvdat[,"pval"] == 1)     # p value = 1 means no subclonal cnv
tmp = sapply(setdiff(id.sub, id.xy), function(i) (2*cnvdat[i,"BAF"]-1)/(cnvdat[i,"nMaj1_A"] - 1 - cnvdat[i,"BAF"]*(cnvdat[i,"nMaj1_A"]+cnvdat[i,"nMin1_A"]-2)))
cellularity0 <- median(tmp[which(abs(tmp) < 1)])

# map mutations to CNV intervals according to their positions
cnv_ssm = vector("list", nrow(cnvdat))
names(cnv_ssm) = rownames(cnvdat)
for(i in 1:length(cnv_ssm)){
  chr = cnvdat[i,"chr"]
  start = cnvdat[i,"startpos"]
  end = cnvdat[i,"endpos"]
  if(chr == "X")
    id <- rownames(ssmdat)[which(ssmdat$CHROM %in% c("X","Y"))]
  else
    id = rownames(ssmdat)[which(ssmdat$CHROM == chr)]
  for(j in id)
    if(ssmdat[j,"POS"] >= start && ssmdat[j,"POS"] <= end)
      cnv_ssm[[i]] = c(cnv_ssm[[i]], j)
}

# get the normal and tumor copy number for each SSM
copyNumber.ssm = matrix(2, nrow = nrow(ssmdat), ncol = 3)
colnames(copyNumber.ssm) = c("ref_cn","var_cn", "coef")
rownames(copyNumber.ssm) = rownames(ssmdat)
# Assuming all samples are from male, i.e., the normal copy number of X nd Y chromosome is 1
id.xy <- c(which(ssmdat$CHROM == "X"),which(ssmdat$CHROM == "Y"))
copyNumber.ssm[id.xy, "ref_cn"] = 1

# Use battenberg CNV calling as the variant copy number
for(i in 1:length(cnv_ssm)){
  cnv1 <- cnvdat[i,"nMaj1_A"] + cnvdat[i,"nMin1_A"]
  frac1 <- cnvdat[i,"frac1_A"]
  if(is.na(cnvdat[i,"frac2_A"]))
    cnv <- cnv1 * frac1
  else{
    cnv2 <- cnvdat[i,"nMaj2_A"] + cnvdat[i,"nMin2_A"]
    frac2 <- cnvdat[i,"frac2_A"]
    cnv <- cnv1*frac1 + cnv2*frac2
  }
  if(cnvdat[i,"chr"] %in% c("X","Y"))
    copyNumber.ssm[cnv_ssm[[i]], "var_cn"] <- cnv/2
  else
    copyNumber.ssm[cnv_ssm[[i]], "var_cn"] = cnv
  # copyNumber.ssm[cnv_ssm[[i]], "var_cn"] = cnv
}

copyNumber.ssm[,"coef"] <- cellularity0*copyNumber.ssm[,"var_cn"] + (1 - cellularity0)*copyNumber.ssm[,"ref_cn"]


# get read counts of mutation reads and total reads
mutReads <- as.numeric(unlist(lapply(strsplit(as.vector(tumor_stat[,"AD"]),','), '[[', 2)))
totalReads <- as.integer(as.vector(tumor_stat[,'DP']))

# Gibbs sampler
GS.data.binomial <- subclone.dirichlet.gibbs(a = mutReads, d = totalReads, copyNumberCoef = copyNumber.ssm[,"coef"])
# density estimator
Gibbs.subclone.density.est(GS.data.binomial)

no.muts <- length(mutReads)
if(no.muts > 1000){
  window.size <- min(50, round(25000/(no.muts-1000)))
  window.size <- max(10, window.size)
} else
  window.size <- 50
  
cluster.assignment <- getClusterAssignment(GS.data.binomial, density.file = "density.txt", window.size = window.size)
#remove empty clusters
occupied.clusters = sort(unique(cluster.assignment$most.likely.cluster), decreasing = TRUE)
no.clusters = length(occupied.clusters)
optima = cluster.assignment$localOptima[occupied.clusters]
optima <- sort(optima, decreasing = TRUE)
assignments = match(cluster.assignment$most.likely.cluster,occupied.clusters)
cellularity = max(optima)
# no.muts = length(assignments)
co.clustering = array(0,c(no.muts,no.muts))
for(c in 1:no.clusters){
  indices = which(assignments==c)
  co.clustering[indices,indices] = 1
}

write.table(cellularity,"GISL_subchallenge1A.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(no.clusters,"GISL_subchallenge1B.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(1:no.clusters,table(assignments),optima),"GISL_subchallenge1C.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(assignments,"GISL_subchallenge2A.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(co.clustering,"GISL_subchallenge2B.txt",row.names=F,col.names=F,quote=F,sep="\t")


######################## Subchallenge 3: Constructing phylogenetic trees ################
nodes <- matrix(nrow = length(optima), ncol = 4)
colnames(nodes) <- c("proportion", "parent","lchild","rchild")
nodes[,"proportion"] <- optima
nodes[1, "parent"] <- 0

unused <- which(is.na(nodes[,"parent"]))
while(length(unused)) {
  leaves <- which(is.na(nodes[,"lchild"]))
  node <- leaves[which.min(nodes[leaves, "parent"])]
  nodes[node,"lchild"] <- unused[1]
  nodes[unused[1], "parent"] <- node
  for(i in unused[-1]) {
    if(nodes[i, "proportion"] < (nodes[node, "proportion"] - nodes[unused[1],"proportion"])){
      nodes[i,"parent"] <- node
      nodes[node, "rchild"] <- i
      break
    }
  }
  unused <- which(is.na(nodes[,"parent"]))
}

# Find ancestors for each cluster
ADM <- matrix(0, nrow = no.muts, ncol = no.muts)
ancestors <- vector("list", no.clusters)
if(no.clusters > 1){
  ancestors[[1]] <- NULL
  for(i in 2:no.clusters){
    parent <- nodes[i, "parent"]
    # parent node's ancestors are also its ancestors
    ancestors[[i]] <- union(parent, ancestors[[parent]])
    id <- which(assignments == i)
    # get ancestors of cluster i
    id.anc <- which(assignments %in% ancestors[[i]])
    ADM[id.anc, id] <- 1
  }
}

write.table(nodes[,"parent"],"GISL_subchallenge3A.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(ADM,"GISL_subchallenge3B.txt",row.names=F,col.names=F,quote=F,sep="\t")
