# In this version, I do prediction
# I split the observations 50:50, select 1 gene randomly to predict, and use the observed values of the remaining genes to do the prediction
library(MegaLMM)
library(BGLR)
library(MCMCglmm)
library(rrBLUP)
library(phenix)
library(lme4qtl)
library(lme4)
source('Estimate_gcor_prediction.R')

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(tidyr)
library(rstan)

seed = as.numeric(commandArgs(t=T)[1])
if(is.na(seed)) seed = 01
seed = seed + 00
n_threads=RcppParallel::defaultNumThreads()-1
registerDoParallel(n_threads)
set.seed(seed)

dataset = 'At_1001'
results_dir = sprintf('Results_%s',dataset)
predictions_dir = sprintf('Predictions_%s',dataset)
try(dir.create(results_dir))
try(dir.create(predictions_dir))

K = read.csv('Data/AT_GE_K.csv',check.names = F,row.names=1)
K = as.matrix(K)
K = K + diag(1e-10,nrow(K))

vst_matrix = fread('Data/At_GE_VST.csv',data.table=F,header = T)
rownames(vst_matrix) = vst_matrix[,1]
vst_matrix = vst_matrix[,-1]
vst_matrix = as.matrix(vst_matrix)
vst_matrix = vst_matrix[,match(rownames(K),colnames(vst_matrix))]
Y_all = t(vst_matrix)

Y_all = scale(Y_all)
# Y_all = Y_all[1:50,]
Y_all = Y_all[sample(1:nrow(Y_all)),]
data_full = data.frame(ID = rownames(Y_all))
K = K[as.character(data_full$ID),as.character(data_full$ID)]

n = nrow(K)
data = data.frame(ID = rownames(K))

tmpdir = sprintf('%s_tmp_%03d_all',dataset,seed)
dir.create(tmpdir)
setwd(tmpdir)

# divide into training, testing
n_test = floor(n/2)
n_train = n-n_test
testing = 1:n_test
training = n_test + 1:n_train

Yn = Y_all[testing,]
Yo = Y_all[training,]
rm(Y_all)

data_train = droplevels(data[training,,drop=F])
data_test = droplevels(data[testing,,drop=F])

Knn = K[testing,testing]
Kno = K[testing,training]
Koo = K[training,training]
sKnn = svd(Knn)
sKoo = svd(Koo)
sK = svd(K)

Koo_inv = solve(Koo)
Kinv_nn_inv = solve(solve(K)[testing,testing])
KnoKoo_inv = Kno %*% Koo_inv
sKinv_nn_inv = svd(Kinv_nn_inv)

# ## MTG2
# # create .fam file
fam = data.frame(FID = data_train$ID,IID = data_train$ID,father=0,mother=0,gender=0,phenotype=0)
write.table(fam,file = 'AT_train.fam',sep=' ',col.names = F,row.names=F,quote=F)

# # write grm file
grm = data.frame(rep(1:n_train,n_train),rep(1:n_train,each = n_train),c(Koo/mean(diag(Koo))))[lower.tri(Koo,diag=T),]
write.table(grm,file = 'AT_train.grm',sep=' ',col.names=F,row.names=F,quote=F)
#
# # do eig
system(sprintf('../Software/mtg2 -p AT_train.fam -g AT_train.grm -pca %d',n_train))

p_test = 1
# ps = 2^c(2:13) + p_test
# ps = 2^c(13) + p_test
ps = c(2^c(5:13)+p_test,ncol(Yn))
# ps = ncol(Yn)
results = c()
g_cor_results = c()
Gs = list()
traits = sample(1:ncol(Yo))
p=4 + p_test
t=1
# traits = c(traits[1],traits)

predictions_mat = c()
Lambda_stats = c()

# r = list()
for(p in ps) {
  predictions = list()
  predictions[['obs']] = Yn[,traits[1]]
  print(p)
  set.seed(seed)
  traits_p = traits[1:p]
  Yp = Yo[,traits_p]
  Yp_joint = rbind(cbind(matrix(NA,nr = n_test,nc = p_test),Yn[,traits_p[-c(1:p_test)]]),Yp)



  if(p < 100) {
    # create phen file
    phen = cbind(fam[,1:2],Yp)
    write.table(phen,file = 'AT_train.phen',sep=' ',col.names = F,row.names=F,quote=F)

    # run MTG2
    time = system.time(mtg_conv <- system(sprintf('../Software/mtg2 -p AT_train.fam -d AT_train.phen -eig AT_train.grm -bv AT_train.bv -out AT_train.out -mod %d -cove 1 -thread %d -nit 1000',p,n_threads),intern = T))
    # fixef = unlist(fread('head -n1 AT_train.bv',data.table=F)[1,]) # On Mac OSX
    fixef = fread('AT_train.bv.fsl',data.table=F,h=T)[,2] # on Unix
    bv = fread('tail -n+2 AT_train.bv',data.table=F,fill=T)#[-1,1:2]
    bv = as.data.frame(tidyr::pivot_wider(cbind(ID = fam[,1],bv),names_from = 'V1',values_from = 'V2'))
    rownames(bv) = bv$ID
    bv = as.matrix(bv[,-1])
    bv = bv[colnames(Koo),,drop=FALSE]

    res = rbind(fread('grep Va AT_train.out',data.table=F),fread('grep cova AT_train.out',data.table=F))
    G_mtg2 = matrix(0,p,p)
    diag(G_mtg2) = res[1:p,2]
    G_mtg2[upper.tri(G_mtg2)] = res[-c(1:p),2]
    G_mtg2[lower.tri(G_mtg2)] = t(G_mtg2)[lower.tri(G_mtg2)]

    res = rbind(fread('grep Ve AT_train.out',data.table=F),fread('grep cove AT_train.out',data.table=F))
    R_mtg2 = matrix(0,p,p)
    diag(R_mtg2) = res[1:p,2]
    R_mtg2[upper.tri(R_mtg2)] = res[-c(1:p),2]
    R_mtg2[lower.tri(R_mtg2)] = t(R_mtg2)[lower.tri(R_mtg2)]

    t = 1
    G12 = G_mtg2[t,-t,drop = F]
    G22 = G_mtg2[-t,-t,drop = F]
    R22 = R_mtg2[-t,-t,drop = F]

    U = sKinv_nn_inv$u
    d = sKinv_nn_inv$d
    UtYhat = t(U) %*% (sweep(Yn[,traits_p][,-t,drop=FALSE],2,fixef[-t],'-') - KnoKoo_inv %*% bv[,-t])
    for(i in 1:n_test) {
      aGR = d[i]*G22 + R22
      UtYhat[i,] = solve(aGR,UtYhat[i,])
    }
    Vci_Yhat = (U) %*% UtYhat
    pred = KnoKoo_inv %*% bv[,t] + Kinv_nn_inv %*% Vci_Yhat %*% t(G12)

    predictions[[paste('MTG2',p-1,sep='.')]] = pred

    results = bind_rows(results,data.frame(Method = 'MTG2',p = p-p_test,time = time[3], converged = sum(grepl('Likelihood not converged >>> may need a longer iterations',mtg_conv)) != 0))
  }

  if(p < 100) {
    p = p
    uty = t(sKoo$u) %*% Yp
    data_train$ut1 = t(sKoo$u) %*% matrix(1,nr = n_train,nc=1)
    D = diag(sKoo$d)
    rownames(D) = colnames(D) = rownames(Koo)
    library(Matrix)
    D = as(D,'dgCMatrix')

    for(i in 1:p) data_train[[paste0('uty',i)]] = uty[,i]
    fixed = formula(paste('cbind(',paste(paste('uty',1:p,sep=''),collapse=','),')~0+ut1:trait'))
    prior = list(R = list(V = diag(1,p),nu = p),G = list(G1 = list(V = diag(1,p),nu = p)))#,alpha.mu = 1,alpha.V = diag(1,p))))
    Dinv = as(diag(1/sKoo$d),'dgCMatrix')
    rownames(Dinv) = colnames(Dinv) = rownames(Koo)
    time = system.time(m4 <- MCMCglmm(fixed,random = ~us(trait):ID,rcov = ~us(trait):units,data = data_train,prior = prior,family = rep('gaussian',p),pl=F,pr=T,
                                      ginverse = list(ID = Dinv),nitt = 7000,burn = 5000,thin = 2))
    m4s = summary(m4)
    G_MCMCglmm = matrix(0,p,p)
    G_MCMCglmm[] = m4s$Gcovariances[,'post.mean']
    R_MCMCglmm = matrix(0,p,p)
    R_MCMCglmm[] = m4s$Rcovariances[,'post.mean']

    U_hat = sKoo$u %*% matrix(colMeans(m4$Sol[,-c(1:p)]),ncol = p)  # re-rotate data
    fixef = colMeans(m4$Sol[,c(1:p)])

    t = 1
    G12 = G_MCMCglmm[t,-t,drop = F]
    G22 = G_MCMCglmm[-t,-t,drop = F]
    R22 = R_MCMCglmm[-t,-t,drop = F]

    U = sKinv_nn_inv$u
    d = sKinv_nn_inv$d
    UtYhat = t(U) %*% (sweep(Yn[,traits_p][,-t,drop=FALSE],2,fixef[-t],'-') - KnoKoo_inv %*% U_hat[,-t])
    for(i in 1:n_test) {
      aGR = d[i]*G22 + R22
      UtYhat[i,] = solve(aGR,UtYhat[i,])
    }
    Vci_Yhat = (U) %*% UtYhat
    pred = KnoKoo_inv %*% U_hat[,t] + Kinv_nn_inv %*% Vci_Yhat %*% t(G12)

    predictions[[paste('MCMCglmm',p-1,sep='.')]] = pred

    ess = apply(m4$Sol[,-c(1:p)],2,ess_bulk)
    results = bind_rows(results,data.frame(Method = 'MCMCglmm',p = p-p_test,time = time[3],matrix(quantile(ess,seq(0,1,by=.1)),nr=1)))
   }

  if(p < n) {
    try({
      time = system.time(m5 <- phenix(Yp_joint,Q = sK$u,lam_K = sK$d))
  
      pred = m5$U[1:n_test,1:p_test,drop=FALSE]
      predictions[[paste('phenix',p-1,sep='.')]] = pred
  
      results = bind_rows(results,data.frame(Method = 'phenix',p = p-p_test,time = time[3]))
    })
  }

  #MegaLMM
  for(k_factor in c(.5,1,2)) {
    k = ceiling(max(4,floor(min(n/4,p/2))) * k_factor)
    run_parameters = MegaLMM_control(
      max_NA_groups = 2,
      save_current_state = FALSE,
      scale_Y = FALSE,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
      h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
      h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
      burn = 0000,  # number of burn in samples before saving posterior samples
      thin = 2,
      K = k # number of factors
    )
    
    priors = MegaLMM_priors(
      tot_Y_var = list(V = 0.5,   nu = 5),      # Prior variance of trait residuals after accounting for fixed effects and factors
      tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
      Lambda_prior = list(
        sampler = sample_Lambda_prec_horseshoe,
        prop_0 = 0.1,
        delta = list(shape = 3, scale = 1),
        delta_iterations_factor = 100
      ),
      h2_priors_resids_fun = function(h2s,n) 1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
      h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
    )
    
    MegaLMM_state = setup_model_MegaLMM(Yp_joint,            # n x p data matrix
                                  ~(1|ID),  # RHS of base model for factors and residuals. Fixed effects defined here only apply to the factor residuals.
                                  data=data_full,         # the data.frame with information for constructing the model matrices
                                  relmat = list(ID = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                                  run_parameters=run_parameters,
                                 run_ID = 'MegaLMM_AT_1'
    )
    maps = make_Missing_data_map(MegaLMM_state)
    MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)

    MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
    MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
    MegaLMM_state = initialize_MegaLMM(MegaLMM_state)
    MegaLMM_state$Posterior$posteriorSample_params = c('Lambda_m_eff')#
    MegaLMM_state$Posterior$posteriorFunctions = list(factor_var = 'get_factor_variance(MegaLMM_state)'
                                                      # ,pred = 'Z[1:n_test,,drop=FALSE] %*% (U_R[,1:p_test,drop=FALSE] + U_F %*% Lambda[,1:p_test,drop=FALSE])'
    )
    MegaLMM_state = clear_Posterior(MegaLMM_state)

    n_samples = 7000
    p2 = min(100,p)
    time = system.time({
      for(i in 1:10) {
        print(i)
        MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
        MegaLMM_state <- sample_MegaLMM(MegaLMM_state,500)
        Lambda_stats = dplyr::bind_rows(Lambda_stats,data.frame(
          trait = colnames(Yo)[traits_p[1]],
          p = p, k_factor = k_factor,k=k,i=i,factor=1:k,
          Lambda_m_eff = MegaLMM::get_posterior_mean(MegaLMM_state$Posterior$Lambda_m_eff)[,1],
          factor_var = MegaLMM::get_posterior_mean(MegaLMM_state$Posterior$factor_var)[,1]
        ))
      }
      MegaLMM_state$Posterior$posteriorSample_params = c('Lambda_m_eff')#
      MegaLMM_state$Posterior$posteriorFunctions = list(factor_var = 'get_factor_variance(MegaLMM_state)'
                                                        ,pred = 'Z[1:n_test,,drop=FALSE] %*% (U_R[,1:p_test,drop=FALSE] + U_F %*% Lambda[,1:p_test,drop=FALSE])'
                                                        ,U_samples = 'U_R[,1:p2] + U_F %*% Lambda[,1:p2]'
      )
      MegaLMM_state = clear_Posterior(MegaLMM_state)
      MegaLMM_state <- sample_MegaLMM(MegaLMM_state,n_samples - MegaLMM_state$current_state$nrun)
      pred_samples = MegaLMM_state$Posterior$pred
      U_samples = MegaLMM_state$Posterior$U_samples
    })
    pred = MegaLMM::get_posterior_mean(pred_samples)
    predictions[[paste(sprintf('MegaLMM_%0.2f_',k_factor),p-1,sep='.')]] = pred
    # system(sprintf('rm -rf %s',MegaLMM_state$Posterior$folder))

    ess = apply(U_samples[,-c(1:n_test),],c(2,3),ess_bulk)

    results = bind_rows(results,data.frame(Method = sprintf('MegaLMM_%0.2f_',k_factor),p = p-p_test,time = time[3],matrix(quantile(c(ess),seq(0,1,by=.1)),nr=1)))
  }

  # rrBLUP
    res_rrBLUP = mixed.solve(Yp_joint[,1],K=K)
    pred = matrix(res_rrBLUP$u[1:n_test],nc=1)

    predictions[[paste('rrBLUP',p-1,sep='.')]] = pred
    
  # collect results
  g_cor_MCMCglmm = foreach(res = predictions[-1],.combine = 'rbind') %dopar% {
    if(!is.null(dim(res))) res = res[,1]
    estimate_gcor(data.frame(ID=data_test$ID,obs = predictions$obs,pred = res),Knn,sKnn,method = 'MCMCglmm',normalize = T)
  }
  
  predictions_methods = sapply(names(predictions),function(x) paste(head(strsplit(x,'.',fixed=T)[[1]],n=-1),collapse='.'))
  predictions_ps = sapply(names(predictions),function(x) tail(strsplit(x,'.',fixed=T)[[1]],n=1))
  new_predictions = which(predictions_ps == p-1)
  new_predictions = na.omit(new_predictions)
  g_cor_results = rbind(g_cor_results,
                      data.frame(method = predictions_methods[new_predictions],
                       p = predictions_ps[new_predictions],
                       trait = colnames(Yo)[traits_p[1]],
                       p_cor = sapply(predictions[new_predictions],function(pred) cor(predictions$obs,c(pred))),
                       g_cor = g_cor_MCMCglmm[,'g_cor'],
                       Rhat = g_cor_MCMCglmm[,'Rhat'])
                )
  write.csv(g_cor_results,file = sprintf('../%s/Gcor_%03d_subset3.csv',results_dir,seed),row.names=F)
  
  write.csv(results,file = sprintf('../%s/Times_%s_subset3.csv',results_dir,tmpdir),row.names=F)
  predictions_mat_new = do.call(cbind,predictions)
  colnames(predictions_mat_new) = names(predictions)
  predictions_mat = cbind(predictions_mat,predictions_mat_new)
  write.csv(predictions_mat,file = sprintf('../%s/predictions_mat_%s_subset3.csv',predictions_dir,tmpdir),row.names = F)
  write.csv(Lambda_stats,file = sprintf('../%s/Lambda_stats_mat_%03d_subset3.csv',results_dir,seed),row.names = F)
}
