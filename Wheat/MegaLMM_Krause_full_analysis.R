library(data.table)
library(rrBLUP)
library(tidyr)
library(MegaLMM)

source('../Method_comparison/Estimate_gcor_prediction.R')

set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)


foldid = 1
foldid = as.numeric(commandArgs(t=T)[1])
if(is.na(foldid)) foldid = 1
set.seed(foldid)


# load data

year = '2014-15'
trt = "Optimal Flat"
results_dir = 'Results_1415_OF_Bgcor'

try(dir.create(results_dir))
geno = fread('dataverse_files/Krause_et_al_2018_Genotypes.csv',data.table=F)
rownames(geno) = geno[,1]
geno = as.matrix(geno[,-1])
# K = A.mat(2*geno-1)
#
# write.csv(K,file = 'K.csv')

K = fread('K.csv',data.table=F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])


BLUEs = fread('dataverse_files/Krause_et_al_2018_Yield_BLUEs.csv',data.table=F)
BLUPs = fread('dataverse_files/Krause_et_al_2018_Yield_iid_BLUPs.csv',data.table=F)


BLUEs = subset(BLUEs,`Breeding Cycle` == year & `Managed_Treatment` == trt)
BLUPs = subset(BLUPs,`Breeding Cycle` == year & `Managed_Treatment` == trt)
BLUEs$GID = as.character(BLUEs$GID)
BLUPs$GID = as.character(BLUPs$GID)
K_year = K[BLUEs$GID,BLUEs$GID]
X = geno[match(BLUEs$GID,rownames(geno)),]
X = X[,colMeans(!is.na(X)) > 0.5 & colMeans(X==0,na.rm=T)>.05]
colMeans_X = colMeans(X,na.rm=T)
X[is.na(X)] = matrix(colMeans_X,nrow = nrow(X),ncol = ncol(X),byrow=T)[is.na(X)]

HTP = fread('dataverse_files/Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv',data.table = F)
HTP = subset(HTP, Breeding_Cycle == year & Managed_Treatment == trt)

HTP_tall = pivot_longer(HTP,cols = Hyper_BLUE_398nm:Hyper_BLUE_847nm,values_to = 'Hyper')
HTP_tall$Phenotype = paste(HTP_tall$name,HTP_tall$Phenotyping_Date,sep='::')

HTP_wide = pivot_wider(HTP_tall,names_from = Phenotype,values_from = Hyper,id_cols = GID)
HTP_wide = HTP_wide[match(BLUEs$GID,HTP_wide$GID),]


# test GP - yes, it's successful. Mean rho ~0.4-0.5. Not sure why they validate with BLUPs instead of BLUEs since they train on BLUEs? But I can do the same.
# The correlation with BLUPs is lower than on BLUEs.

data = data.frame(BLUE = BLUEs$Grain_Yield_BLUE,
                  BLUP = BLUPs$Grain_Yield_iid_BLUP,
                  GID = BLUEs$GID,
                  GID2 = BLUEs$GID)

print(foldid)
nas = sample(1:nrow(data),.5*nrow(data))
data$yNA = data$BLUE
data$yNA[nas] = NA

runID = sprintf('MegaLMM_Krause_K_%d',foldid)


  # run MegaLMM
  run_parameters = MegaLMM_control(
    # num_NA_groups = 0,
    drop0_tol = 1e-10,
    scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
    simulation = FALSE, # Are you running against simulated data (ex from a call to new_halfSib_simulation above)? If so, you can provide the setup list and it will make some QC plots
    h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
    h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
    burn = 00,  # number of burn in samples before saving posterior samples
    K = 100 # number of factors
  )

  priors = MegaLMM_priors(
    tot_Y_var = list(V = 0.5,   nu = 10),      # Prior variance of trait residuals after accounting for fixed effects and factors
    tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
    Lambda_prior = list(                     # Prior fo rfactor loadings
      sampler = sample_Lambda_prec_ARD,
      Lambda_df = 3,
      delta_1 = list(shape = 20, rate = 1/2),
      delta_2 = list(shape = 3, rate = 1)
    ),
    h2_priors_resids_fun = function(h2s,n)  1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
    h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
  )

  # data matrix with testing data masked
  Y = cbind(data$yNA,HTP_wide[,-1])

  # scaling is not necesary, but makes visualization easier
  Y = scale(Y)

  MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                ~ 1 + (1|GID),   # in this case, only genotype ID used as a predictor
                                data=data,         # the data.frame with information for constructing the model matrices
                                relmat = list(GID = K_year), # covariance matrices for the random effects. If not provided, assume uncorrelated
                                run_parameters=run_parameters,
                                run_ID = runID
  )
  # column_groups = unname(sapply(colnames(Y_BLUP),function(x) strsplit(x,'::')[[1]][2]))
  maps = make_Missing_data_map(MegaLMM_state,2,verbose=T)    # depending on structure of missing data, might increase the 2nd argument
  MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)
  MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)  # get starting values

  MegaLMM_state = initialize_MegaLMM(MegaLMM_state)   # pre-calculate some matrix inverses
  MegaLMM_state = clear_Posterior(MegaLMM_state)

  MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','U_F','U_R','Eta','F_h2')#,'Eta_mean')   # these are the ones to save for standard models
  # MegaLMM_state$Posterior$posteriorSample_params = c('Eta')#,'Eta_mean')
  # MegaLMM_state$Posterior$posteriorMean_params = c()
  MegaLMM_state = clear_Posterior(MegaLMM_state)   # start with a clean (empty) set of posterior samples


  # for the MCMC, I generally break the chain into a set of chunks, so the model can be saved, checked, and adjusted throughout, rather than waiting until the end
  n_iter = 100;  # how many samples to collect at once?
  system(sprintf('rm %s/U_pred_samples.csv',runID))
  for(i  in 1:40) {
    print(sprintf('Run %d',i))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter,grainSize=1)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)

    # assuming the focal trait is in column 1, this extracts the predictions of the masked data from just this trait and saves it to a file.
    # These predictions are not saved by default
    U = get_posterior_FUN(MegaLMM_state,U_R[,1:2,drop=F] + U_F %*% t(Lambda[1:2,,drop=F]))
    U = U[,nas,1]
    fwrite(as.data.table(U),file = sprintf('%s/U_pred_samples.csv',runID),append=T)

    MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
    print(MegaLMM_state) # print status of current chain
    plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf

    # set of commands to run during burn-in period to help chain converge
    if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i <= 20) {
      MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
      MegaLMM_state = clear_Posterior(MegaLMM_state)
      print(MegaLMM_state$run_parameters$burn)
    }
  }

  # now that the chain is done, we can re-load the posterior samples, and summarize
  MegaLMM_state$Posterior = readRDS(sprintf('%s/Posterior/Posterior_base.rds',runID))
  # U = get_posterior_mean(MegaLMM_state,U_R + U_F %*% t(Lambda),bychunk = T)
  U = get_posterior_mean(MegaLMM_state,U_R[,1:2,drop=F] + U_F %*% t(Lambda[1:2,,drop=F]),bychunk = T)


  # code for estimating prediction accuracy. First is pearson correlation which is biased if traits are measured on the same plots. The second is the more robust method, but takes much longer.
  results = data.frame(Method = 'MegaLMM',
                       pearson = cor(data$BLUP[nas],U[nas,1])/sqrt(h2_BLUE),
                       g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = U[nas,1]),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']])
  results$fold = foldid
  write.csv(results,file = sprintf('%s/results_MegaLMM_fold_%d.csv',results_dir,foldid))


