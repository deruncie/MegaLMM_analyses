library(data.table)
library(DESeq2)

exp = fread('Data/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv',data.table=F,check.names = F)

rownames(exp) = exp[,1]
exp = exp[,-1]
ACC_ID = colnames(exp)
ACC_ID = sub('X','',ACC_ID)
colnames(exp) = ACC_ID
exp = as.matrix(exp)
exp = exp[rowMeans(exp)>=10,]

# write.table(cbind(ACC_ID,ACC_ID),file = 'ACC_ID.txt',row.names=F,quote=F,col.names = F)

# load kinship matrix from 1001 genomes
# download file from: http://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5
K_1001_accessions = fread('Data/1001_accessions.csv')
K_1001 = fread('Data/1001_kinship.csv')
K_1001 = as.matrix(K_1001[,-1])
rownames(K_1001) = colnames(K_1001) = K_1001_accessions$x

exp = exp[,colnames(exp) %in% K_1001_accessions$x]

K = K_1001[colnames(exp),colnames(exp)]

normalize_K = function(K){
  # centers K and then normalizes to have mean diagonal = 1
  n = nrow(K)
  J = matrix(1,n,1)
  M = diag(1,n) - J %*% solve(t(J) %*% J) %*% t(J)
  centered_K = M %*% K %*% t(M)
  centered_K = centered_K/mean(diag(centered_K))
  rownames(centered_K) = rownames(K)
  colnames(centered_K) = colnames(K)
  centered_K
}
K = normalize_K(K)
write.csv(K,file = 'Data/AT_GE_K.csv',row.names=T,col.names = T)

# normalize data

colData = data.frame(ACC_ID = colnames(exp))
design = model.matrix(~1,colData)

dds = DESeqDataSetFromMatrix(countData = exp,colData = colData,design=~1)
dds_vst = varianceStabilizingTransformation(dds,blind=T)
vst_matrix = assay(dds_vst)
fwrite(as.data.frame(vst_matrix),file = 'Data/At_GE_VST.csv',row.names = T)
