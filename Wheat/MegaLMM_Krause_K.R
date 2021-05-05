library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4qtl)
library(lme4)
library(MegaLMM)
source('Estimate_gcor_prediction.R')
MegaLMM::set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)
results_dir = 'Results'

trial = as.numeric(commandArgs(t=T)[1])
if(is.na(trial)) trial = 1
foldid = as.numeric(commandArgs(t=T)[2])
if(is.na(foldid)) foldid = 1
k = as.numeric(commandArgs(t=T)[3])
if(is.na(k)) k = 100

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

Knn = K[nas,nas]
sKnn = svd(Knn)

## rrBLUP - K
res_K = mixed.solve(BLUEs$yNA,K = K)
h2_BLUE = res_K$Vu/(res_K$Vu+res_K$Ve)

# runID = sprintf('/group/runciegrp/Projects/MegaLMM/Krause/MegaLMM_Krause_K_%d',foldid)
runID = sprintf('MegaLMM_Krause_K_%d',foldid)


# predict_MegaLMM = function(data,HTP_wide,K_year,runID,nas) {
  # run MegaLMM
  run_parameters = MegaLMM_control(
    drop0_tol = 1e-10,
    scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
    h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
    h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
    burn = 00,  # number of burn in samples before saving posterior samples
    K = k # number of factors
  )
  
  priors = MegaLMM_priors(
    tot_Y_var = list(V = 0.5,   nu = 10),      # Prior variance of trait residuals after accounting for fixed effects and factors
    tot_F_var = list(V = 18/20, nu = 100000),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
    Lambda_prior = list(
      sampler = sample_Lambda_prec_horseshoe,
      prop_0 = 0.1,
      delta = list(shape = 3, scale = 1),
      delta_iterations_factor = 100
    ),
    h2_priors_resids_fun = function(h2s,n) 1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
    h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
  )
  
  Y = cbind(GY = BLUEs$yNA,HTP)

  Y = scale(Y)

  MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                      ~ 1 + (1|GID),
                                      data=BLUEs,         # the data.frame with information for constructing the model matrices
                                      relmat = list(GID = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                                      run_parameters=run_parameters,
                                      run_ID = runID
  )
  maps = make_Missing_data_map(MegaLMM_state,2,verbose=T)
  MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)
  MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
  MegaLMM_state = initialize_MegaLMM(MegaLMM_state)
  
  # MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','U_F','U_R','Eta','F_h2')#,'Eta_mean')
  MegaLMM_state$Posterior$posteriorSample_params = c('Lambda')
  MegaLMM_state$Posterior$posteriorFunctions = list(pred = 'U_R[,1] + U_F %*% Lambda[,1]')
  MegaLMM_state = clear_Posterior(MegaLMM_state)

  n_iter = 100;  # how many samples to collect at once?
  for(i  in 1:40) {
    print(sprintf('Run %d',i))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
    
    MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
    print(MegaLMM_state) # print status of current chain
    plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
    
    # set of commands to run during burn-in period to help chain converge
    if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i <= 20) {
      MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
      print(MegaLMM_state$run_parameters$burn)
    }
  }
  
  U = get_posterior_mean(load_posterior_param(MegaLMM_state,'pred'))

  results = data.frame(Method = 'MegaLMM_K',
                       pearson = cor(BLUEs$BLUP[nas],U[nas,1])/sqrt(h2_BLUE),
                       g_cor = estimate_gcor(data.frame(ID=BLUEs$GID[nas],obs = BLUEs$BLUP[nas],pred = U[nas,1]),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']])
  results$fold = foldid
  write.csv(results,file = sprintf('%s/results_MegaLMM_fold_%d.csv',results_dir,foldid))


