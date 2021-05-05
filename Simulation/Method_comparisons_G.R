# create simulated data with known G and R (based on correlations in gene expression dataset)
# assess accuracy of each method at estimating G and R
library(MegaLMM)
library(MCMCglmm)
library(phenix)

library(data.table)
library(dplyr)
library(tidyr)
library(doParallel)

seed = as.numeric(commandArgs(t=T)[1])
if(is.na(seed)) seed = 101
seed = seed + 00
n_threads=RcppParallel::defaultNumThreads()-1
registerDoParallel(n_threads)
set.seed(seed)

dataset = 'At_1001'
results_dir = sprintf('Results_G_accuracy_%s',dataset)
try(dir.create(results_dir))

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
Y_all = Y_all[sample(1:nrow(Y_all)),]
data_full = data.frame(ID = rownames(Y_all))
K = K[as.character(data_full$ID),as.character(data_full$ID)]
sK = svd(K)

n = nrow(K)
data = data.frame(ID = rownames(K))

tmpdir = sprintf('%s_G_tmp_%03d_all',dataset,seed)
dir.create(tmpdir)
setwd(tmpdir)

# ## MTG2
# # create .fam file
fam = data.frame(FID = data$ID,IID = data$ID,father=0,mother=0,gender=0,phenotype=0)
write.table(fam,file = 'AT_train.fam',sep=' ',col.names = F,row.names=F,quote=F)

# # write grm file
grm = data.frame(rep(1:n,n),rep(1:n,each = n),c(K/mean(diag(K))))[lower.tri(K,diag=T),]
write.table(grm,file = 'AT_train.grm',sep=' ',col.names=F,row.names=F,quote=F)
#
# # do eig
system(sprintf('../Software/mtg2 -p AT_train.fam -g AT_train.grm -pca %d',n))


# create random G and R matrices
p = 128
G = cor(Y_all[,sample(1:ncol(Y_all),p)])
R = cor(Y_all[,sample(1:ncol(Y_all),p)])

h2s = runif(p,.1,.8)
G = sqrt(h2s) * G %*% diag(sqrt(h2s))
R= sqrt(1-h2s) * R %*% diag(sqrt(1-h2s))

# create data
U = t(chol(K)) %*% rstdnorm_mat(nrow(K),ncol(G)) %*% chol(G)
E = rstdnorm_mat(nrow(K),ncol(G)) %*% (chol(R))

Y = U+E

get_results = function(G,G_hat,name){
  tri = upper.tri(G,diag=F)
  results = data.frame(RMSE_var = sqrt(mean((diag(G_hat)-diag(G))^2)),
             RMSE_G = sqrt(mean((G_hat[tri] - G[tri])^2)),
             cor_covs = cor(c(G_hat[tri]),c(G[tri])))
  colnames(results) = paste(name,colnames(results),sep='_')
  results
}

# estimate R for subsets
ps = 2^c(2:7)

results = c()

for(p in ps) {
  print(p)
  Yp = Y[,1:p]
  
  if(p < 100) {
    # MTG2
    # run MTG2
    phen = cbind(fam[,1:2],Yp)
    write.table(phen,file = 'AT_train.phen',sep=' ',col.names = F,row.names=F,quote=F)
    
    mtg_conv <- system(sprintf('../Software/mtg2 -p AT_train.fam -d AT_train.phen -eig AT_train.grm -bv AT_train.bv -out AT_train.out -mod %d -cove 1 -thread %d -nit 1000',p,n_threads))
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
    
    results = bind_rows(results,data.frame(
      Method = 'MTG2',p = p,seed=seed,
      get_results(G[1:p,1:p],G_mtg2,'G'),
      get_results(R[1:p,1:p],R_mtg2,'R')
    ))
    
    # MCMCglmm
    library(MCMCglmm)
    p2 = p
    uty = t(sK$u) %*% Yp
    data$ut1 = t(sK$u) %*% matrix(1,nr = n,nc=1)
    D = diag(sK$d)
    rownames(D) = colnames(D) = rownames(K)
    library(Matrix)
    D = as(D,'dgCMatrix')
    
    for(i in 1:p) data[[paste0('uty',i)]] = uty[,i]
    # fixed = formula(paste('cbind(',paste(paste('y',1:p2,sep=''),collapse=','),')~0+trait'))
    fixed = formula(paste('cbind(',paste(paste('uty',1:p2,sep=''),collapse=','),')~0+ut1'))
    prior = list(R = list(V = diag(1,p2),nu = p2),G = list(G1 = list(V = diag(1,p2),nu = p2)))#,alpha.mu = 1,alpha.V = diag(1,p2))))
    # Ainv = as(solve(A),'dgCMatrix')
    Kinv = as(diag(1/sK$d),'dgCMatrix')
    rownames(Kinv) = colnames(Kinv) = rownames(K)
    m4 <- MCMCglmm(fixed,random = ~us(trait):ID,rcov = ~us(trait):units,data = data,prior = prior,family = rep('gaussian',p2),ginverse = list(ID = Kinv),nitt = 7000,burn = 5000,thin = 1)
    m4s = summary(m4)
    G_MCMCglmm = matrix(0,p2,p2)
    G_MCMCglmm[] = m4s$Gcovariances[,'post.mean']
    G_MCMCglmm_ESS = G_MCMCglmm
    G_MCMCglmm_ESS[] = m4s$Gcovariances[,'eff.samp']
    R_MCMCglmm = matrix(0,p2,p2)
    R_MCMCglmm[] = m4s$Rcovariances[,'post.mean']
    R_MCMCglmm_ESS = R_MCMCglmm
    R_MCMCglmm_ESS[] = m4s$Rcovariances[,'eff.samp']
    
    results = bind_rows(results,data.frame(
      Method = 'MCMCglmm',p = p,seed=seed,
      get_results(G[1:p,1:p],G_MCMCglmm,'G'),
      get_results(R[1:p,1:p],R_MCMCglmm,'R')
    ))
  }
  
  # phenix
    m5 <- phenix(Yp,Q = sK$u,lam_K = sK$d)
    G_phenix = m5$B
    R_phenix = m5$E
    
    results = bind_rows(results,data.frame(
      Method = 'phenix',p = p,seed=seed,
      get_results(G[1:p,1:p],G_phenix,'G'),
      get_results(R[1:p,1:p],R_phenix,'R')
    ))
  
  # MegaLMM
  
    k = min(64,2*p)
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
    
    MegaLMM_state = setup_model_MegaLMM(Yp,            # n x p data matrix
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
    MegaLMM_state$Posterior$posteriorSample_params = c()
    MegaLMM_state$Posterior$posteriorMean_params = c()
    
    n_samples = 7000
    for(i in 1:10) {
      print(i)
      MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 1) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
      MegaLMM_state <- sample_MegaLMM(MegaLMM_state,500)
    }
    # MegaLMM_state$Posterior$posteriorSample_params = c('Lambda_m_eff')#
    MegaLMM_state$Posterior$posteriorFunctions = list(G = 't(Lambda) %*% ((F_h2[1,])*Lambda) + diag(resid_h2[1,]/tot_Eta_prec[1,])',
                                                      R = 't(Lambda) %*% ((1-F_h2[1,])*Lambda) + diag((1-resid_h2[1,])/tot_Eta_prec[1,])'
    )
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    MegaLMM_state <- sample_MegaLMM(MegaLMM_state,n_samples - MegaLMM_state$current_state$nrun)
    G_MegaLMM = MegaLMM::get_posterior_mean(MegaLMM_state$Posterior$G)
    R_MegaLMM = MegaLMM::get_posterior_mean(MegaLMM_state$Posterior$R)
    
    results = bind_rows(results,data.frame(
      Method = 'MegaLMM',p = p,seed=seed,
      get_results(G[1:p,1:p],G_MegaLMM,'G'),
      get_results(R[1:p,1:p],R_MegaLMM,'R')
    ))
  
  write.csv(results,file = sprintf('../%s/G_accuracy_%03d.csv',results_dir,seed))
}


# library(cowplot)
# library(ggplot2)
# library(dplyr)
# library(plyr)
# library(tidyr)
# 
# dodge_width = .5
# (p1 <- ggplot(results,aes(x=log2(p),y=G_RMSE_G)) + ylab('RMSE of genetic covariances') + xlab('# traits') + #ggtitle('G') + 
#     # geom_point(aes(color = Method,group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
#     geom_pointrange(fatten=2,stat="summary", fun.data="mean_se",aes(color = Method,group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
#     geom_smooth(aes(color=Method,group=Method),se=F) +
#     scale_x_continuous(breaks = unique(log2(results$p)),labels = unique(results$p)) +
#     theme_cowplot() + background_grid(major = 'xy')
# )
# 
# (p2 <- ggplot(results,aes(x=log2(p),y=R_RMSE_G)) + ylab('RMSE of residual covariances') + xlab('# traits')  + #ggtitle('R') + 
#     geom_pointrange(fatten=2,stat="summary", fun.data="mean_se",aes(color = Method,group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
#     geom_smooth(aes(color=Method,group=Method),se=F) +
#     scale_x_continuous(breaks = unique(log2(results$p)),labels = unique(results$p)) +
#     theme_cowplot() + background_grid(major = 'xy')
# )
# leg = get_legend(p1 + theme(legend.position = 'bottom',legend.justification = 'center'))
# nl = theme(legend.position = 'none')
# save_plot(plot_grid(plot_grid(p1+nl,p2+nl,labels = c('A','B')),leg,nrow = 2,rel_heights = c(1,.2)),file = 'G_accuracy.pdf',base_height = 5)
