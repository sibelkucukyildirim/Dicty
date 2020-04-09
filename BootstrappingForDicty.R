#Created by Megan Behringer
library(boot)

Dicty_Mutations<-read.table("/usr/local//DictyTableS1.txt", sep= "\t", header=TRUE)
#Make value for total mutations per line
Dicty_Mutations$TotalMut<-Dicty_Mutations$GC_to_AT+Dicty_Mutations$AT_to_GC+Dicty_Mutations$AT_to_TA+Dicty_Mutations$GC_to_TA+Dicty_Mutations$AT_to_CG+Dicty_Mutations$GC_to_CG
#Make value for total number of GC sites queried per line
Dicty_Mutations$GCSites<-Dicty_Mutations$Num_C_sites+Dicty_Mutations$Num_G_sites
#Make value for total number of AT sites queried per line
Dicty_Mutations$ATSites<-Dicty_Mutations$Num_A_sites+Dicty_Mutations$Num_T_sites
#Make value for total indels
Dicty_Mutations$TotalIndel<-Dicty_Mutations$Insertions+Dicty_Mutations$Deletions

#Split data set into Lynch and Queller subets
Dicty_Lynch <- Dicty_Mutations[which(Dicty_Mutations$Lab=='Lynch'),]
Dicty_QS <- Dicty_Mutations[which(Dicty_Mutations$Lab=='Q_S'),]

#Make functions to determine different mutation rates (statistics) for bootstrapping d=dataset w=weighted values; values are weighted because of variance in number of sites and number of generations
MutationRate <-function(d,w) sum(d$TotalMut*w)/(sum(d$Gens*w)*sum(d$Total_Sites*w))
GC_to_ATRate <-function(d,w) sum(d$GC_to_AT*w)/(sum(d$Gens*w)*sum(d$GCSites*w))
AT_to_GCRate <-function(d,w) sum(d$AT_to_GC*w)/(sum(d$Gens*w)*sum(d$ATSites*w))
AT_to_TARate <-function(d,w) sum(d$AT_to_TA*w)/(sum(d$Gens*w)*sum(d$ATSites*w))
GC_to_TARate <-function(d,w) sum(d$GC_to_TA*w)/(sum(d$Gens*w)*sum(d$GCSites*w))
AT_to_CGRate <-function(d,w) sum(d$AT_to_CG*w)/(sum(d$Gens*w)*sum(d$ATSites*w))
GC_to_CGRate <-function(d,w) sum(d$GC_to_AT*w)/(sum(d$Gens*w)*sum(d$GCSites*w))
InsertionRate<-function(d,w) sum(d$Insertions*w)/(sum(d$Gens*w)*sum(d$Total_Sites*w))
DeletionRate<-function(d,w) sum(d$Deletions*w)/(sum(d$Gens*w)*sum(d$Total_Sites*w))
IndelRate<-function(d,w) sum(d$TotalIndel*w)/(sum(d$Gens*w)*sum(d$Total_Sites*w))

#command to run single group boostraping boot(d,statistic, R=<Number of Bootstraps, stpye="<weighting variable from statistic>")
#Make and edit these commants to get the specific mutation rates above
MutationRateAll.boot<-boot(Dicty_Mutations, MutationRate, R = 1000, stype = "w")
MutationRateLynch.boot<-boot(Dicty_Lynch, MutationRate, R = 1000, stype = "w")
MutationRateQS.boot<-boot(Dicty_QS, MutationRate, R = 1000, stype = "w")

#command to get confidence intervals using BCa method
#Make and edit these commands to get CIs from each mutation rate from each lab.
coords<-cbind(rbind(c(1,4),c(1,5)))
colnames(coords)<-c("low","high")
MutationRateAll.CI<-boot.ci(MutationRateAll.boot, conf = 0.95,type =  "bca")
#command to only print the last two columns of the boot.ci out put above (low.CI, high.CI)
MutationRateAll.CI$bca[coords]

