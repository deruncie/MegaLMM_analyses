library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4qtl)
library(lme4)

foldid = as.numeric(commandArgs(t=T)[1])
if(is.na(foldid)) foldid = 1
set.seed(foldid)

source('data_prep_Krause.R')
source('../Method_comparison/Estimate_gcor_prediction.R')


# res = mixed.solve(data$yNA,K=K_year)
# m1 = relmatLmer(yNA~(1|GID),data,relmat = list(GID = K_year))
Knn = K_year[nas,nas]
Kno = K_year[nas,-nas]
Koo = K_year[-nas,-nas]
sKnn = svd(Knn)

Hnn = K_HTP[nas,nas]
Hno = K_HTP[nas,-nas]
Hoo = K_HTP[-nas,-nas]

Hbest_nn = K6_HTP[nas,nas]
Hbest_no = K6_HTP[nas,-nas]
Hbest_oo = K6_HTP[-nas,-nas]

# u = as.matrix(t(m1@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m1)[[1]])[rownames(Koo),1])
# pred_lme4qtl_G = Kno %*% solve(Koo,u)

sKoo = svd(Koo)
Q = sKoo$u
D = diag(sKoo$d)
rownames(D) = colnames(D) = rownames(Koo)
QtX = t(Q) %*% model.matrix(~1,data[-nas,])
m1 = relmatLmer(t(Q) %*% yNA~0+QtX+(1|GID),subset(data,!is.na(yNA)),relmat = list(GID = D))
u = Q %*% as.matrix(t(m1@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m1)[[1]])[rownames(Koo),1])
pred_lme4qtl_G = Kno %*% solve(Koo,u)

sHoo = svd(Hoo)
Q = sHoo$u
D = diag(sHoo$d)
rownames(D) = colnames(D) = rownames(Hoo)
QtX = t(Q) %*% model.matrix(~1,data[-nas,])
m1 = relmatLmer(t(Q) %*% yNA~0+QtX+(1|GID),subset(data,!is.na(yNA)),relmat = list(GID = D))
u = Q %*% as.matrix(t(m1@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m1)[[1]])[rownames(Hoo),1])
pred_lme4qtl_H = Hno %*% solve(Hoo,u)

sHbest_oo = svd(Hbest_oo)
Q = sHbest_oo$u
D = diag(sHbest_oo$d)
rownames(D) = colnames(D) = rownames(Hbest_oo)
QtX = t(Q) %*% model.matrix(~1,data[-nas,])
m1 = relmatLmer(t(Q) %*% yNA~0+QtX+(1|GID),subset(data,!is.na(yNA)),relmat = list(GID = D))
u = Q %*% as.matrix(t(m1@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m1)[[1]])[rownames(Hbest_oo),1])
pred_lme4qtl_Hbest = Hbest_no %*% MASS::ginv(Hbest_oo) %*% u



# verifying by direct calculation
# nn = nrow(data)-length(nas)
# vars = as.data.frame(VarCorr(m1))$vcov
# plot(pred_lme4qtl_G,vars[1]*Kno %*% solve(vars[1]*Koo+vars[2]*diag(1,nn),subset(data,!is.na(yNA))$yNA - fixef(m1)));abline(0,1)

m2 = relmatLmer(yNA~(1|GID)+(1|GID2),data[-nas,],relmat = list(GID = K_year,GID2 = K_HTP)) #K_HTP
u1 = as.matrix(t(m2@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m2)[[1]]))[rownames(Koo),1]
u2 = as.matrix(t(m2@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(m2)[[2]]))[rownames(Hoo),1]
pred_lme4qtl_KH = Kno %*% solve(Koo,u1) + Hno %*% MASS::ginv(Hoo) %*% u2

# verifying by direct calculation
# vars = as.data.frame(VarCorr(m2))$vcov
# Vc = vars[1]*Koo+vars[2]*Hoo + vars[3]*diag(1,nn)
# plot(pred_lme4qtl_KH-(vars[1]*Kno+vars[2]*Hno) %*% solve(Vc,subset(data,!is.na(yNA))$yNA - fixef(m2)));abline(0,1)


m3 = relmatLmer(yNA~(1|GID)+(1|GID2),data[-nas,],relmat = list(GID = K_year,GID2 = K6_HTP))
Hnn = K6_HTP[nas,nas]
Hno = K6_HTP[nas,-nas]
Hoo = K6_HTP[-nas,-nas]
u1 = as.matrix(t(m3@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(m3)[[1]]))[rownames(Koo),1]
u2 = as.matrix(t(m3@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(m3)[[2]]))[rownames(Hoo),1]
pred_lme4qtl_KH6 = Kno %*% solve(Koo,u1) + Hno %*% MASS::ginv(Hoo) %*% u2


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
             pearson = cor(data$BLUP[nas],pred_lme4qtl_G)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = pred_lme4qtl_G),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'GBLUP_H',
             pearson = cor(data$BLUP[nas],pred_lme4qtl_H)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = pred_lme4qtl_H),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'GBLUP_Hbest',
             pearson = cor(data$BLUP[nas],pred_lme4qtl_Hbest)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = pred_lme4qtl_Hbest),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'GBLUP_KH',
             pearson = cor(data$BLUP[nas],pred_lme4qtl_KH)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = pred_lme4qtl_KH),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]),
  data.frame(Method = 'GBLUP_KHbest',
             pearson = cor(data$BLUP[nas],pred_lme4qtl_KH6)/sqrt(h2_BLUE),
             g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = pred_lme4qtl_KH6),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']])
)

results$fold = foldid
write.csv(results,file = sprintf('%s/results_H_fold_%d.csv',results_dir,foldid))
