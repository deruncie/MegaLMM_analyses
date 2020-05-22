library(data.table)
library(lme4qtl)
library(lme4)

# uses lme4 to get BLUPs for each variety across all trials, with some variety:trial combos masked.
# then uses lme4qtl to fit genetic values for each variety


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

YNA_tall = reshape2::melt(YNA)
YNA_tall = na.omit(YNA_tall)
m1 = lmer(value~Var2 + (1|Var1),YNA_tall)
blups = as.data.frame(ranef(m1)$Var1)
K = A_mat[rownames(blups),rownames(blups)]
sK = svd(K)
blups$qtX = c(t(sK$u) %*% model.matrix(~1,blups))
blups$qty = t(sK$u) %*% blups[,1]
blups$ID = rownames(blups)
D = diag(sK$d)
rownames(D) = colnames(D) = rownames(K)
res = relmatLmer(qty~0+qtX+(1|ID),blups,relmat = list(ID = D))
u_hat = sK$u %*% (as.matrix(t(res@optinfo$relmat$relfac$ID) %*% ranef(res)$ID[,1])[blups$ID,1])

U = matrix(NA,nrow = nrow(Y),ncol = ncol(Y),dimnames = list(rownames(Y),colnames(Y)))
U[blups$ID,] = u_hat

alternatives_results = list(
  Y = Y,
  mask = mask,
  U_mean = U
)
try(dir.create('G2F_byTrait_results'))
saveRDS(alternatives_results,file = sprintf('G2F_byTrait_results/mean_%s_%d.rds',trait,seed))
