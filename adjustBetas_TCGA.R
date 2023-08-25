#####======================================================================#####
### Correct TCGA BRCA top5000 beta values for infiltration
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncoloy and Pathology

##Git for data download and preprocessing of TCGA BRCA data (Staaf-lab pipeline)
#https://github.com/StaafLab/tcgaBrca
#resluting data object can be obtained on demand from MA/Staaf-lab to avoid extensive download/run times

##Git for beta adjustment
#https://github.com/StaafLab/adjustBetas

##Reduced data set for demo purposes

################################################################################
##Set directory paths

##set/create own home directory below:
HOME<-"~/Documents/Duke/lab/adjustBetas"

################################################################################
##load required packages

if(!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel") }

library(doParallel)

if(!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel") }

library(parallel)

if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer") }

library(RColorBrewer)

if(!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap") }

library(pheatmap)

##source - flexmix loaded on source
source(paste0(HOME,"/function_correctBetas.r"))

################################################################################
### Load reduced data set, do correction and plot heatmap

##data set reduced to top 5000 CpGs in set of 630 TCGA BRCA tumors
#load("workspace_tcgaBrca_top5000.RData")

ls()
# [1] "adjustBeta"        "betaUnadjusted"    "cpgAnnotations"    "purityVector"     
# [5] "sampleAnnotations"

##Function for adjusting beta values
#adjustBeta

##Top 5000 CpG probes
str(betaUnadjusted)
 # num [1:5000, 1:630] 0.341 0.057 0.012 0.022 0 0.004 0.021 0.892 0.022 0.163 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:5000] "cg09248054" "cg25340711" "cg06443533" "cg16601494" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

##custom CpG annotations
#cpgAnnotations

##Sample purities from TCGA
str(purityVector)
 # Named num [1:630] 0.85 0.73 0.69 0.62 0.64 0.53 0.89 0.47 0.56 0.65 ...
 # - attr(*, "names")= chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

##TCGA sample annotations
#sampleAnnotations

################################################################################
###correct betas using multicore 

#need to feed seed to function so fully reproducible

no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-1,"cores","\n")

cl <- makeCluster(no_cores-1)  
registerDoParallel(cl)  

##estimated runtime is ~500 sec on 7 cores

clusterEvalQ(cl, {
  library("flexmix")
})


indir <- '/Users/danielmonyak/Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/Monyak Data/CpG Site Selection'
#proj_dir <- file.path(indir, 'adjustBetas')
#outdir <- proj_dir
proj_dir <- file.path(indir, 'desmedt')
outdir <- file.path(indir, 'desmedt', 'adjustBetas')

#beta_values_balanced_CpGs = read.table(file.path(indir, 'beta_values_ALL_balanced_CpGs_NONANS.txt'), sep='\t', header=TRUE, row.names=1)
beta_values_balanced_CpGs = read.table(file.path(proj_dir, 'adjustBetas', 'beta_values_balanced_CpGs_NONANS.txt'), sep='\t', header=TRUE, row.names=1)

#CPE_purity = read.table('/Users/danielmonyak/Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/Monyak Data/CpG Site Selection/CPE_purity.txt', sep='\t', row.names=1, header=TRUE)
CPE_purity = read.table(file.path(proj_dir, 'LUMP_purity.txt'), sep='\t', row.names=1, header=TRUE)
good.samples = rownames(CPE_purity)[!is.na(CPE_purity)]
CPE_purity_vector = CPE_purity[good.samples, ]
names(CPE_purity_vector) = make.names(good.samples)

#purity <- purity_vector
#data <- betaUnadjusted[1:100, ]
purity <- CPE_purity_vector
data <- beta_values_balanced_CpGs[, names(purity)]

##add rng seed to each row and pass to each paralell instance
betaRun<-cbind(seed=1:nrow(data),data)
betaNames<-colnames(data)

#clusterSetRNGStream(cl, 20200918) ##using this way to pass seed will not make exactly replicable..
res<-parRapply(cl = cl, betaRun, adjustBeta,purity=purity,snames=betaNames,seed=TRUE)

#res <- adjustBeta(methylation=data[1, ],purity=purity,snames=betaNames,nmax=3,nrep=3,seed=FALSE)

res2<-parRapply(cl = cl, betaRun, adjustBeta,purity=purity,snames=betaNames,seed=TRUE)

table(unlist(lapply(res,function(x) x$n.groups)))
  #  2    3 
  # 20 4980 

table(unlist(lapply(res,function(x) x$n.groups)),unlist(lapply(res2,function(x) x$n.groups)))
  #      2    3
  # 2   20    0
  # 3    0 4980

table( unlist(lapply(1:length(res),function(x) { all( unlist(res[[x]],use.names=FALSE) == unlist(res[[x]],use.names=FALSE) ) }) ) )
# TRUE 
# 5000 
table( unlist(lapply(1:length(res),function(x) { all( unlist(res[[x]],use.names=FALSE) == unlist(res2[[x]],use.names=FALSE) ) }) ) )
# TRUE 
# 5000 

##same results across runs.
rm(betaRun,betaNames,res2,cl,no_cores)

################################################################################
###plot top5k clusters - adjuster order

##check stats
beta_values_TUMOR_balanced_CpGs <- do.call("rbind",lapply(res,function(x) x$y.tum))
beta_values_NORMAL_balanced_CpGs <- do.call("rbind",lapply(res,function(x) x$y.norm))
#temp3<-do.call("rbind",lapply(res,function(x) x$y.orig))


write.table(beta_values_TUMOR_balanced_CpGs, file.path(outdir, 'beta_values_TUMOR_balanced_CpGs.txt'), quote=FALSE, sep='\t')
write.table(beta_values_NORMAL_balanced_CpGs, file.path(outdir, 'beta_values_NORMAL_balanced_CpGs.txt'), quote=FALSE, sep='\t')
