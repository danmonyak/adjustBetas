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
load("workspace_tcgaBrca_top5000.RData")

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

#TCGA_beta_values_balanced_CpGs = read.table('/Users/danielmonyak/Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/Monyak Data/CpG Site Selection/ALL/beta_values_ALL_balanced_CpGs.txt', sep='\t', header=TRUE, row.names=1)
TCGA_beta_values_balanced_CpGs = read.table('/Users/danielmonyak/Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/Monyak Data/CpG Site Selection/beta_values_ALL_balanced_CpGs_NONANS.txt', sep='\t', header=TRUE, row.names=1)

CPE_purity = read.table('/Users/danielmonyak/Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/Monyak Data/CpG Site Selection/CPE_purity.txt', sep='\t', row.names=1, header=TRUE)
good.samples = rownames(CPE_purity)[!is.na(CPE_purity)]
CPE_purity_vector = CPE_purity[good.samples, ]
names(CPE_purity_vector) = make.names(good.samples)

#purity <- purity_vector
#data <- betaUnadjusted[1:100, ]
purity <- CPE_purity_vector
data <- TCGA_beta_values_balanced_CpGs[, names(purity)]

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
temp1<-do.call("rbind",lapply(res,function(x) x$y.tum))
temp2<-do.call("rbind",lapply(res,function(x) x$y.norm))
temp3<-do.call("rbind",lapply(res,function(x) x$y.orig))

all(abs(temp3 - data) < 0.001)
                              
table(apply(temp1,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp2,1,function(x) sum(is.na(x))))
#   0
#5000
table(apply(temp3,1,function(x) sum(is.na(x))))
#   0
#5000

#####
sample.stdevs.orig <- apply(data, 2, sd)
sample.stdevs.tum <- apply(temp1, 2, sd)
sample.stdevs.norm <- apply(temp2, 2, sd)

stdevs.df <- as.data.frame(cbind(orig=sample.stdevs.orig, tum=sample.stdevs.tum, norm=sample.stdevs.norm))
            
sample.means.orig <- apply(data, 2, mean)
sample.means.tum <- apply(temp1, 2, mean)
sample.means.norm <- apply(temp2, 2, mean)

means.df <- as.data.frame(cbind(orig=sample.means.orig, tum=sample.means.tum, norm=sample.means.norm))
 
ggplot(stdevs.df, aes(x = orig, y = tum)) + geom_point(size=1)
ggplot(stdevs.df, aes(x = orig, y = norm)) + geom_point(size=1)

ggplot(means.df, aes(x = orig, y = tum)) + geom_point(size=1)
hist(means.df$norm)
#####

            
table(cpgAnnotations[rownames(temp1),"chr"])
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
#   545   269   212   254   155   147    99   170   273    68   205   362   129 
# chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
#    42    54   252   191   371   398   389   357    58 

##do clustering rows+columns
c1<-cutree( hclust( as.dist( 1-cor(temp1) ),method="ward.D"),5)
c2<-hclust( as.dist( 1-cor(temp1) ),method="ward.D")
c3<-cutree( hclust( as.dist( 1-cor(temp3) ),method="ward.D"),5)
r1<-hclust( dist(temp1),method="ward.D")

sample_anno<-data.frame(adj5000=as.character(c1),
  unadj5000=as.character(c3),
  ER=sampleAnnotations$ER,
  PR=sampleAnnotations$PR,
  HER2=sampleAnnotations$HER2,
  TNBC=as.character(as.integer(sampleAnnotations$TNBC)),
  PAM50=sampleAnnotations$pam50.full,stringsAsFactors=FALSE
  )
rownames(sample_anno)<-colnames(temp1)
sample_anno<-sample_anno[,ncol(sample_anno):1]

my_colour = list(unadj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    adj5000=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A","4"="#984EA3","5"="#FF7F00"),
    ER = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    PR = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A"
      ),
    HER2 = c("[Not Available]"="#FFFF33",
      "[Not Evaluated]"="#FF7F00",
      "Equivocal"="#377EB8",
      "Indeterminate"="#984EA3",
      "Negative"="#E41A1C",
      "Positive"="#4DAF4A",
      "NA"="white"
      ),
    TNBC = c("1"="black","0"="lightgrey"),
    PAM50 = c("Basal"="#E41A1C","LumB"="#377EB8","LumA"="#4DAF4A","Her2"="#984EA3","Normal"="#808080","NA"="white"),stringsAsFactors=FALSE
  )

tiff("top5k_heatmap_pear_eucl_adjClust_adjBeta.tiff",width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp1,cluster_rows = r1, cluster_cols = c2
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, adj data, adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff("top5k_heatmap_pear_eucl_adjClust_normBeta.tiff",width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp2,cluster_rows = r1, cluster_cols = c2
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, adj data, adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

tiff("top5k_heatmap_pear_eucl_adjClust_unadjBeta.tiff",width=10*500,height=13*500,units="px",res=500,compression="lzw")
pheatmap(temp3,cluster_rows = r1, cluster_cols = c2
  ,show_rownames=F,show_colnames=F
  ,main="top 5000 by sd, unadj data, adj clust , pearC/euclR",cutree_cols=5
  ,annotation_col=sample_anno,annotation_colors=my_colour
)
dev.off()

##END
