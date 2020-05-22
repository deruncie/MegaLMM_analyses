library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4qtl)
library(lme4)
library(MegaLMM)
library(MCMCglmm)
# source('MegaLMM_Krause.R')
# source('MegaLMM_Krause_geno.R')
# set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)


foldid = as.numeric(commandArgs(t=T)[1])
if(is.na(foldid)) foldid = 1
set.seed(foldid)

year = '2014-15'
trt = "Optimal Flat"
results_dir = 'Results_1415_OF_Bgcor'

# year = '2013-14'
# trt = "Severe Drought"
# results_dir = 'Results_1314_SD'

try(dir.create(results_dir))



BLUEs = fread('Data/Krause_et_al_2018_Yield_BLUEs.csv',data.table=F)
BLUPs = fread('Data/Krause_et_al_2018_Yield_iid_BLUPs.csv',data.table=F)


BLUEs = subset(BLUEs,`Breeding Cycle` == year & `Managed_Treatment` == trt)
BLUPs = subset(BLUPs,`Breeding Cycle` == year & `Managed_Treatment` == trt)
BLUEs$GID = as.character(BLUEs$GID)
BLUPs$GID = as.character(BLUPs$GID)

geno = fread('Data/Krause_et_al_2018_Genotypes.csv',data.table=F)
rownames(geno) = geno[,1]
geno = as.matrix(geno[,-1])

# drop individuals with < 100 markers
# then markers with MAF < 0.05 or <30% genotyped markers
# geno = geno[rowSums(!is.na(geno)) >= 100,]
# BLUEs = subset(BLUEs,GID %in% rownames(geno))
# BLUPs = subset(BLUPs,GID %in% rownames(geno))

geno = geno[BLUEs$GID,]

# K = A.mat(2*geno-1)
# write.csv(K,file = 'Data/K.csv')
# # 
# # restrict markers to those with <50% missing data, impute remaining missing data before calculating distance matrix
# geno = geno[,colMeans(!is.na(geno)) >  0.5]
# geno[is.na(geno)] = matrix(colMeans(geno,na.rm=T),nr = nrow(geno),nc = ncol(geno),byrow=T)[is.na(geno)]
# dim(geno)
# write.csv(geno,'Data/geno_year.csv')
# D = as.matrix(dist(geno))
# write.csv(D,file = 'Data/D.csv')

K = fread('Data/K.csv',data.table=F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

D = fread('Data/D.csv',data.table=F,h=T)
rownames(D) = D[,1]
D = as.matrix(D[,-1])

X = fread('Data/geno_year.csv',data.table=F,h=T)
rownames(X) = X[,1]
X = as.matrix(X[,-1])


# HTP contains 10 timepoints with 60 wavelengths per timepoint. HTP1 groups this into 4 growth stages
HTP = fread('Data/Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv',data.table = F)
HTP = subset(HTP, Breeding_Cycle == year & Managed_Treatment == trt)

HTP_tall = pivot_longer(HTP,cols = Hyper_BLUE_398nm:Hyper_BLUE_847nm,values_to = 'Hyper')
HTP_tall$Phenotype = paste(HTP_tall$name,HTP_tall$Phenotyping_Date,sep='::')

HTP_wide = pivot_wider(HTP_tall,names_from = Phenotype,values_from = Hyper,id_cols = GID)
HTP_wide = HTP_wide[match(BLUEs$GID,HTP_wide$GID),]

K_HTP = tcrossprod(scale(HTP_wide[,-1])) / ncol(HTP_wide[,-1])
rownames(K_HTP) = colnames(K_HTP) = HTP_wide$GID

data = data.frame(BLUE = BLUEs$Grain_Yield_BLUE,
                  BLUP = BLUPs$Grain_Yield_iid_BLUP,
                  GID = BLUEs$GID,
                  GID2 = BLUEs$GID)



  print(foldid)
  nas = sample(1:nrow(data),.5*nrow(data))
  data$yNA = data$BLUE
  data$yNA[nas] = NA
 
  Knn = K[nas,nas]
  Kno = K[nas,-nas]
  Koo = K[-nas,-nas]
  sKnn = svd(Knn)

 
  sKoo = svd(Koo)
  Q = sKoo$u
  D = diag(sKoo$d)
  rownames(D) = colnames(D) = rownames(Koo)
  QtX = t(Q) %*% model.matrix(~1,data[-nas,])
  m1 = relmatLmer(t(Q) %*% yNA~0+QtX+(1|GID),subset(data,!is.na(yNA)),relmat = list(GID = D))
  u = Q %*% as.matrix(t(m1@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m1)[[1]])[rownames(Koo),1])
  pred_lme4qtl_G = Kno %*% solve(Koo,u)

  vars1 = as.data.frame(VarCorr(m1))
  h2_BLUE = vars1$vcov[1]/sum(vars1$vcov)

