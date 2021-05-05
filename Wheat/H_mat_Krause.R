library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4qtl)
library(lme4)
source('Estimate_gcor_prediction.R')
results_dir = 'Results'

trial = as.numeric(commandArgs(t=T)[1])
if(is.na(trial)) trial = 1
foldid = as.numeric(commandArgs(t=T)[2])
if(is.na(foldid)) foldid = 1

folder = sprintf('Trial_%02d',trial)
setwd(folder)

# load data
BLUEs = fread('BLUES.csv',data.table=F)

K = fread('K.csv',data.table=F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

D = fread('D.csv',data.table=F,h=T)
rownames(D) = D[,1]
D = as.matrix(D[,-1])

X = fread('X.csv',data.table=F)
rownames(X) = X[,1]
X = as.matrix(X[,-1])

HTP = fread('HTP.csv',data.table=F)
rownames(HTP) = HTP[,1]
HTP = as.matrix(HTP[,-1])

# select training
set.seed(foldid)
BLUEs$partition = sample(c('Training','Validation'),size=nrow(BLUEs),replace=T)
BLUEs$yNA = BLUEs$Grain_Yield_BLUE
nas = BLUEs$partition == 'Validation'
BLUEs$yNA[nas] = NA
BLUEs$GID2 = BLUEs$GID

H = tcrossprod(scale(HTP)) / ncol(HTP)
rownames(H) = colnames(H) = rownames(HTP)


# res = mixed.solve(data$yNA,K=K_year)
# m1 = relmatLmer(yNA~(1|GID),BLUEs,relmat = list(GID = K))
Knn = K[nas,nas]
Kno = K[nas,!nas]
Koo = K[!nas,!nas]
sKnn = svd(Knn)

Hnn = H[nas,nas]
Hno = H[nas,!nas]
Hoo = H[!nas,!nas]


# u = as.matrix(t(m1@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m1)[[1]])[rownames(Koo),1])
# pred_lme4qtl_G = Kno %*% solve(Koo,u)

sKoo = svd(Koo)
Q = sKoo$u
D = diag(sKoo$d)
rownames(D) = colnames(D) = rownames(Koo)
QtX = t(Q) %*% model.matrix(~1,BLUEs[!nas,])
m1 = relmatLmer(t(Q) %*% yNA~0+QtX+(1|GID),subset(BLUEs,!is.na(yNA)),relmat = list(GID = D))
u = Q %*% as.matrix(t(m1@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m1)[[1]])[rownames(Koo),1])
pred_lme4qtl_G = Kno %*% MASS::ginv(Koo) %*% u

sHoo = svd(Hoo)
Q = sHoo$u
D = diag(sHoo$d)
rownames(D) = colnames(D) = rownames(Hoo)
QtX = t(Q) %*% model.matrix(~1,BLUEs[!nas,])
m1 = relmatLmer(t(Q) %*% yNA~0+QtX+(1|GID),subset(BLUEs,!is.na(yNA)),relmat = list(GID = D))
u = Q %*% as.matrix(t(m1@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m1)[[1]])[rownames(Hoo),1])
# pred_lme4qtl_H = Hno %*% solve(Hoo,u)
pred_lme4qtl_H = Hno %*% MASS::ginv(Hoo) %*% u


# verifying by direct calculation
# nn = nrow(BLUEs)-sum(nas)
# vars = as.data.frame(VarCorr(m1))$vcov
# plot(pred_lme4qtl_G,vars[1]*Kno %*% solve(vars[1]*Koo+vars[2]*diag(1,nn),subset(BLUEs,!is.na(yNA))$yNA - fixef(m1)));abline(0,1)
# plot(pred_lme4qtl_H,vars[1]*Hno %*% solve(vars[1]*Hoo+vars[2]*diag(1,nn),subset(BLUEs,!is.na(yNA))$yNA - fixef(m1)));abline(0,1)

m2 = relmatLmer(yNA~(1|GID)+(1|GID2),BLUEs[!nas,],relmat = list(GID = K,GID2 = H)) #K_HTP
u1 = as.matrix(t(m2@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m2)[[1]]))[rownames(Koo),1]
u2 = as.matrix(t(m2@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(m2)[[2]]))[rownames(Hoo),1]
pred_lme4qtl_KH = Kno %*% MASS::ginv(Koo) %*% u1 + Hno %*% MASS::ginv(Hoo) %*% u2

# verifying by direct calculation
# vars = as.data.frame(VarCorr(m2))$vcov
# Vc = vars[1]*Koo+vars[2]*Hoo + vars[3]*diag(1,nn)
# plot(pred_lme4qtl_KH-(vars[1]*Kno+vars[2]*Hno) %*% solve(Vc,subset(BLUEs,!is.na(yNA))$yNA - fixef(m2)))

vars1 = as.data.frame(VarCorr(m1))
h2_BLUE = vars1$vcov[1]/sum(vars1$vcov)


# results = rbind(
#   data.frame(Method = 'GBLUP',
#              pearson = cor(data$BLUP[nas],pred_lme4qtl_G)/sqrt(h2_BLUE),
#              g_cor = calc_gcor(data,nas,pred_lme4qtl_G,sKnn)),
#   data.frame(Method = 'GBLUP_H',
#              pearson = cor(data$BLUP[nas],pred_lme4qtl_H)/sqrt(h2_BLUE),
#              g_cor = calc_gcor(data,nas,pred_lme4qtl_H,sKnn)),
#   data.frame(Method = 'GBLUP_Hbest',
#              pearson = cor(data$BLUP[nas],pred_lme4qtl_Hbest)/sqrt(h2_BLUE),
#              g_cor = calc_gcor(data,nas,pred_lme4qtl_Hbest,sKnn)),
#   data.frame(Method = 'GBLUP_KH',
#              pearson = cor(data$BLUP[nas],pred_lme4qtl_KH)/sqrt(h2_BLUE),
#              g_cor = calc_gcor(data,nas,pred_lme4qtl_KH,sKnn)),
#   data.frame(Method = 'GBLUP_KHbest',
#              pearson = cor(data$BLUP[nas],pred_lme4qtl_KH6)/sqrt(h2_BLUE),
#              g_cor = calc_gcor(data,nas,pred_lme4qtl_KH6,sKnn))
# )

results = rbind(
  data.frame(Method = 'GBLUP',
             pearson = cor(BLUEs$BLUP[nas],pred_lme4qtl_G)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=BLUEs$GID[nas],obs = BLUEs$BLUP[nas],pred = pred_lme4qtl_G),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'GBLUP_H',
             pearson = cor(BLUEs$BLUP[nas],pred_lme4qtl_H)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=BLUEs$GID[nas],obs = BLUEs$BLUP[nas],pred = pred_lme4qtl_H),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'GBLUP_KH',
             pearson = cor(BLUEs$BLUP[nas],pred_lme4qtl_KH)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=BLUEs$GID[nas],obs = BLUEs$BLUP[nas],pred = pred_lme4qtl_KH),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']])
)

results$fold = foldid
write.csv(results,file = sprintf('%s/results_H_fold_%d.csv',results_dir,foldid))
