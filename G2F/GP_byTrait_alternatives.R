library(phenix)
library(rrBLUP)
library(data.table)
library(foreach)
library(doParallel)

# runs rrBLUP for univariate predictions and phenix for multivariate predictions

id = as.numeric(commandArgs(t=T)[1])
if(is.na(id)) id = 103

dataset = 0+id %/% 100
seed = id %% 100

set.seed(seed)

files = list.files(path='BLUP_matrices/',pattern = 'csv')
Y = read.csv(sprintf('BLUP_matrices/%s',files[dataset]),row.names=1)
data = data.frame(Pedigree_taxa = rownames(Y),stringsAsFactors = F)

trait = sub('_BLUPs.csv','',files[dataset])
trait = gsub(' ','_',trimws(trait))


A_mat = fread('Data/g2f_2014_zeaGBSv27_CenteredIBS_allYears.txt',data.table=F)
rownames(A_mat) = A_mat[,1]
A_mat = as.matrix(A_mat[,-1])
colnames(A_mat) = rownames(A_mat)
A_mat = A_mat/mean(diag(A_mat))

all(rownames(Y) %in% rownames(A_mat))

mask = matrix(F,nrow = nrow(Y),ncol = ncol(Y))
for(i in 1:ncol(Y)) {
  obs = which(!is.na(Y[,i]))
  mask[sample(obs,.2*length(obs)),i] = T
}

YNA = Y
YNA[mask] = NA
YNA = as.matrix(YNA)

A_mat = A_mat[data$Pedigree_taxa,data$Pedigree_taxa]
registerDoParallel(RcppParallel::defaultNumThreads()-1)
U_rrBLUP = foreach(i=1:ncol(Y),.combine = 'cbind') %dopar% {
  res = mixed.solve(YNA[,i],K=A_mat)
  res$u
}

A_mat = A_mat+diag(1e-6,nrow(A_mat))
res_phenix = phenix(YNA,K = A_mat)


alternatives_results = list(
  Y = Y,
  mask = mask,
  U_rrBLUP = U_rrBLUP,
  U_phenix = res_phenix$U,
  Y_phenix = res_phenix$imp
)
try(dir.create('G2F_byTrait_results'))
saveRDS(alternatives_results,file = sprintf('G2F_byTrait_results/alternatives_%s_%d.rds',trait,seed))
