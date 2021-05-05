library(MegaLMM)
library(data.table)
set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)

id = as.numeric(commandArgs(t=T)[1])
if(is.na(id)) id = 101

dataset = 0+id %/% 100
seed = id %% 100

files = list.files(path='BLUP_matrices/',pattern = 'csv')
Y = read.csv(sprintf('BLUP_matrices/%s',files[dataset]),row.names=1)
Y = Y[,apply(Y,2,var,na.rm=T)>0.05]
data = data.frame(Pedigree_taxa = rownames(Y),stringsAsFactors = F)

trait = sub('_BLUPs.csv','',files[dataset])
trait = gsub(' ','_',trimws(trait))


A_mat = fread('Data/g2f_2014_zeaGBSv27_CenteredIBS_allYears.txt',data.table=F)
rownames(A_mat) = A_mat[,1]
A_mat = as.matrix(A_mat[,-1])
colnames(A_mat) = rownames(A_mat)

all(rownames(Y) %in% rownames(A_mat))

runID = sprintf('/group/runciegrp/Projects/MegaLMM/G2F/MegaLMM_full_%s_%d',dataset,seed)
# runID = sprintf('MegaLMM_%s_%d',dataset,seed)

# run MegaLMM
run_parameters = MegaLMM_control(
  # num_NA_groups = 0,
  drop0_tol = 1e-10,
  scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
  simulation = FALSE, # Are you running against simulated data (ex from a call to new_halfSib_simulation above)? If so, you can provide the setup list and it will make some QC plots
  h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = 00,  # number of burn in samples before saving posterior samples
  thin = 2,
  K = 50 # number of factors
)

priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5,   nu = 10),      # Prior variance of trait residuals after accounting for fixed effects and factors
  tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
  Lambda_prior = list(
    sampler = sample_Lambda_prec_horseshoe,
    prop_0 = 0.1,
    delta = list(shape = 3, scale = 1),
    delta_iterations_factor = 100
  ),
  # Lambda_prior = list(
  #   sampler = sample_Lambda_prec_ARD,
  #   Lambda_df = 3,
  #   delta_1 = list(shape = 20, rate = 1/2),
  #   delta_2 = list(shape = 3, rate = 1)
  # ),
  # B2_prior = list(
  #   sampler = sample_B2_prec_horseshoe,
  #   prop_0 = 0.1
  # ),
  # cis_effects_prior = list(
  #   prec = 1
  # ),
  h2_priors_resids_fun = function(h2s,n)  1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
)

MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                              ~ 1 + (1|Pedigree_taxa),
                              data=data,         # the data.frame with information for constructing the model matrices
                              relmat = list(Pedigree_taxa = A_mat), # covariance matrices for the random effects. If not provided, assume uncorrelated
                              run_parameters=run_parameters,
                              run_ID = runID
)
maps = make_Missing_data_map(MegaLMM_state,max_NA_groups = 20,verbose=T)
MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map_list[[4]])
MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state_base = MegaLMM_state
saveRDS(MegaLMM_state_base,file = sprintf('%s/MegaLMM_state_base.rds',runID))
MegaLMM_state = initialize_MegaLMM(MegaLMM_state)
MegaLMM_state = clear_Posterior(MegaLMM_state)

MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2','tot_Eta_prec','F','U_F')
MegaLMM_state$Posterior$posteriorMean_params = c('Eta','Eta_mean')
MegaLMM_state = clear_Posterior(MegaLMM_state)


n_iter = 100;  # how many samples to collect at once?
system(sprintf('rm %s/U_pred_samples.csv',runID))
for(i  in 1:70) {
  print(sprintf('Run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)

  MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
  print(MegaLMM_state) # print status of current chain
  plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf

  # set of commands to run during burn-in period to help chain converge
  if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i < 50) {
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    print(MegaLMM_state$run_parameters$burn)
  }
}

# id = 801
# dataset = 0+id %/% 100
# seed = id %% 100
# files = list.files(path='BLUP_matrices/',pattern = 'csv')
# trait = sub('_BLUPs.csv','',files[dataset])
# trait = gsub(' ','_',trimws(trait))

MegaLMM_state = readRDS(sprintf('/group/runciegrp/Projects/MegaLMM/G2F/MegaLMM_full_%s_%d/MegaLMM_state_base.rds',dataset,seed))
MegaLMM_state$Posterior = readRDS(sprintf('%s/Posterior/Posterior_base.rds',MegaLMM_state$run_ID))
MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)
Gcor_samples = get_posterior_FUN(MegaLMM_state,cov2cor(t(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])))
Gcor = get_posterior_mean(Gcor_samples)
try(dir.create('G2F_byTrait_Gcors'))
saveRDS(Gcor,file = sprintf('G2F_byTrait_Gcors/Gcor_%s_%d.rds',trait,seed))

h2s = get_posterior_mean(MegaLMM_state,{
  gvar = colSums(F_h2[1,]*t(Lambda)^2 + resid_h2[1,]/tot_Eta_prec[1,])
  evar = colSums((1-F_h2[1,])*t(Lambda)^2 + (1-resid_h2[1,])/tot_Eta_prec[1,])
  gvar/(gvar+evar)
})
saveRDS(h2s,file = sprintf('G2F_byTrait_Gcors/h2s_%s_%d.rds',trait,seed))

Lambda = get_posterior_mean(MegaLMM_state,Lambda,bychunk = T)
saveRDS(Lambda,file = sprintf('G2F_byTrait_Gcors/Lambda_%s_%d.rds',trait,seed))

F = get_posterior_mean(MegaLMM_state,F,bychunk = T)
saveRDS(F,file = sprintf('G2F_byTrait_Gcors/F_%s_%d.rds',trait,seed))
U_F = get_posterior_mean(MegaLMM_state,U_F,bychunk = T)
saveRDS(U_F,file = sprintf('G2F_byTrait_Gcors/U_F_%s_%d.rds',trait,seed))

