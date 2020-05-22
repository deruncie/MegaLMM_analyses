library(data.table)
library(lme4)
library(tidyr)


MIN_GENO_PER_TRIAL = 50
MIN_VAR_BLUP = 0.01
nfolds = 20
mask_perc = .2
set.seed(1)

geno_info = fread('Data/g2f_2014_zeaGBSv27_DNAtoTaxa.txt',data.table=F)
pheno = fread('Data/g2f_2014_2017_hybrid_data_clean.csv',data.table=F)
pheno = subset(pheno,!is.na(Tester))
pheno$RepBlock = paste(pheno$Replicate,pheno$Block,sep='.')
pheno$Range = factor(pheno$Range)
pheno$Pass = factor(pheno$Pass)

pheno = subset(pheno,Tester %in% geno_info$DNAName & DNAName_mod %in% geno_info$DNAName)

pheno$Tester_taxa = geno_info$Taxa[match(pheno$Tester,geno_info$DNAName)]
pheno$Male_taxa = geno_info$Taxa[match(pheno$DNAName_mod,geno_info$DNAName)]
pheno$Pedigree_taxa = paste(pheno$Male_taxa,pheno$Tester_taxa,sep='/')
pheno$SiteYear = paste(pheno$`Field-Location`,pheno$Year,sep='-')
# colnames(pheno)[43] = 'ASI [days]'

hybrid_info = unique(pheno[,c('Male_taxa','Tester_taxa','Pedigree_taxa')])

write.table(hybrid_info[,c('Male_taxa','Tester_taxa')],file = 'Data/Hybrid_genotypes.txt',sep='\t',row.names=F,quote=F)

A_mat = fread('Data/g2f_2014_zeaGBSv27_CenteredIBS_allYears.txt',data.table = F)
rownames(A_mat) = A_mat[,1]
A_mat = as.matrix(A_mat[,-1])
colnames(A_mat) = rownames(A_mat)
sum(is.na(A_mat))

sapply(2014:2017,function(year) {
  mean(unique(subset(pheno,Year==year)$Pedigree_taxa) %in% rownames(A_mat))
})

trait_cols = colnames(pheno)[c(26:29,33:36,43)]
SYs = unique(pheno$SiteYear)

# make univariate BLUP matrix
trait. = trait_cols[1]
try(dir.create('BLUP_matrices'))
for(trait. in trait_cols[c(2,3,8,9)]) {
  Y_BLUP_trait = matrix(NA,nr = nrow(hybrid_info),nc = length(SYs))
  rownames(Y_BLUP_trait) = hybrid_info$Pedigree_taxa
  colnames(Y_BLUP_trait) = SYs
  sy = SYs[1]
  for(sy in SYs) {
    sy_data = subset(pheno,SiteYear == sy)
    sy_data$y = sy_data[[trait.]]
    sy_data$y = scale(sy_data$y)
    sy_data = subset(sy_data,!is.na(y))
    if(nrow(sy_data) == 0) next
    # print(sprintf('%s %02f',sy,mean(table(sy_data$Pedigree_taxa)>1)))
    if(mean(table(sy_data$Pedigree_taxa)>1) < 0.25) {
      means = aggregate(y~Pedigree_taxa,sy_data,FUN=mean)
      Y_BLUP_trait[means$Pedigree_taxa,sy] = means$V1
    } else{
        res = lmer(y~RepBlock + (1|Pedigree_taxa),sy_data)
        blups = as.matrix(ranef(res)$Pedigree_taxa)
        Y_BLUP_trait[rownames(blups),sy] = blups[,1]
    }
  }
  Y_BLUP_trait = Y_BLUP_trait[,colSums(!is.na(Y_BLUP_trait)) >= MIN_GENO_PER_TRIAL]
  Y_BLUP_trait = Y_BLUP_trait[,apply(Y_BLUP_trait,2,var,na.rm=T) > MIN_VAR_BLUP]
  print(dim(Y_BLUP_trait))
  write.csv(Y_BLUP_trait,file = sprintf('BLUP_matrices/%s_BLUPs.csv',strsplit(trait.,' [',fixed=T)[[1]][1]),row.names=T)
}



for(file in list.files(path='BLUP_matrices',pattern = 'csv',full.names = T)){
  r = fread(file,data.table=F)
  print(file)
  print(dim(r))
  print(mean(colMeans(!is.na(r[,-1]))))
}

