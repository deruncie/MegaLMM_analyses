library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4qtl)
library(lme4)
library(MegaLMM)
set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)
foldid = 1
# foldid = as.numeric(commandArgs(t=T)[1])
# if(is.na(foldid)) foldid = 1
set.seed(foldid)


source('data_prep_Krause.R')
foldid = as.numeric(commandArgs(t=T)[1])
if(is.na(foldid)) foldid = 1
set.seed(foldid)

# runID = sprintf('/group/runciegrp/Projects/MegaLMM/Krause/MegaLMM_Krause_K_%d',foldid)
runID = sprintf('MegaLMM_Krause_K_full_%d',foldid)


# predict_MegaLMM = function(data,HTP_wide,K_year,runID,nas) {
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
    Lambda_prior = list(
      sampler = sample_Lambda_prec_horseshoe,
      prop_0 = 0.1,
      delta = list(shape = 3, scale = 1),
      delta_iterations_factor = 100
    ),
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

  Y = cbind(data$BLUE,HTP_wide[,-1])
  # Y = cbind(data$yNA,HTP_wide[,-1])

  Y = scale(Y)
  # Y = Y[,1:3]

  MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                ~ 1 + (1|GID),
                                data=data,         # the data.frame with information for constructing the model matrices
                                relmat = list(GID = K_year), # covariance matrices for the random effects. If not provided, assume uncorrelated
                                run_parameters=run_parameters,
                                run_ID = runID
  )
  # column_groups = unname(sapply(colnames(Y_BLUP),function(x) strsplit(x,'::')[[1]][2]))
  maps = make_Missing_data_map(MegaLMM_state,2,verbose=T)
  MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)
  MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
  MegaLMM_state_base = MegaLMM_state
  saveRDS(MegaLMM_state_base,file = sprintf('%s/MegaLMM_state_base.rds',runID))
  MegaLMM_state = initialize_MegaLMM(MegaLMM_state)
  MegaLMM_state = clear_Posterior(MegaLMM_state)

  MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','U_F','U_R','Eta','F_h2','resid_h2')#,'Eta_mean')
  # MegaLMM_state$Posterior$posteriorSample_params = c('Eta')#,'Eta_mean')
  # MegaLMM_state$Posterior$posteriorMean_params = c()
  MegaLMM_state = clear_Posterior(MegaLMM_state)


  n_iter = 100;  # how many samples to collect at once?
  # system(sprintf('rm %s/U_pred_samples.csv',runID))
  for(i  in 1:40) {
    print(sprintf('Run %d',i))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter,grainSize=1)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)

    # U = get_posterior_FUN(MegaLMM_state,U_R[,1:2,drop=F] + U_F %*% t(Lambda[1:2,,drop=F]))
    # U = U[,nas,1]
    # fwrite(as.data.table(U),file = sprintf('%s/U_pred_samples.csv',runID),append=T)

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

  MegaLMM_state$Posterior$Lambda = load_posterior_param(MegaLMM_state,'Lambda')
  MegaLMM_state$Posterior$F_h2 = load_posterior_param(MegaLMM_state,'F_h2')
  MegaLMM_state$Posterior$resid_h2 = load_posterior_param(MegaLMM_state,'resid_h2')
  MegaLMM_state$Posterior$tot_Eta_prec = load_posterior_param(MegaLMM_state,'tot_Eta_prec')
  G_cov_samples = get_posterior_FUN(MegaLMM_state,Lambda %*% diag(F_h2[1,]) %*% Lambda[1,])
  diag_G_samples = get_posterior_FUN(MegaLMM_state,colSums(F_h2[1,]*t(Lambda^2)) + resid_h2[1,]/tot_Eta_prec)
  G_cor_samples = G_cov_samples[,,1] / sqrt(diag_G_samples[,1,] * diag_G_samples[,1,1])
  saveRDS(G_cor_samples,file = sprintf('%s/G_cor_samples.rds',runID))

  G_cor_samples = readRDS('MegaLMM_Krause_K_full_1/G_cor_samples.rds')
  G_cor_mean = colMeans(G_cor_samples)
  G_cor_HPD = get_posterior_HPDinterval(G_cor_samples)
  P_cor = cor(data$BLUE,HTP_wide[,-1])[1,]
  results = data.frame(trait = names(P_cor),G_cor_low = G_cor_HPD[1,-1],G_cor_high = G_cor_HPD[2,-1],G_cor_mean = G_cor_mean[-1],P_cor = P_cor)
  results = separate(results,'trait',c('wavelength','date'),sep='::')
  results$wavelength = as.numeric(substr(results$wavelength,12,14))
  results$date = factor(results$date,levels = unique(results$date),labels = format(as.Date(unique(results$date),format = '%y%m%d'),'%d-%b'))
  (p <- ggplot(results,aes(x=wavelength)) +
    facet_wrap(~date) +
    geom_hline(yintercept = 0,size=.25) +
    geom_ribbon(aes(ymin = G_cor_low,ymax = G_cor_high),alpha = 0.3) +
    geom_line(aes(y = G_cor_mean,col = 'G')) +
    geom_line(aes(y = P_cor,col = 'P')) +
    scale_color_manual(values=c('G' = 'red','P' = 'black'),name = 'Correlation') +
    xlab('Wavelength') + ylab('Correlation') +
    theme(legend.position = c(.8,.15)))
  save_plot('MegaLMM_Krause_K_full_1/wavelength_cors.pdf',p)


  # G2 = get_posterior_FUN(MegaLMM_state,cov2cor(Lambda[1:10,] %*% diag(F_h2[1,]) %*% t(Lambda[1:10,]) + diag(resid_h2[1,1:10]/tot_Eta_prec[1,1:10])))


  # rm(MegaLMM_state)
  # gc()
  #
  # # MegaLMM_state = list()
  # MegaLMM_state = MegaLMM_state_base
  # MegaLMM_state$Posterior = readRDS(sprintf('%s/Posterior/Posterior_base.rds',runID))
  # # U = get_posterior_mean(MegaLMM_state,U_R + U_F %*% t(Lambda),bychunk = T)
  # U = get_posterior_mean(MegaLMM_state,U_R[,1:2,drop=F] + U_F %*% t(Lambda[1:2,,drop=F]),bychunk = T)

  # return(U[nas,1])
  # Eta_samples = load_posterior_param(MegaLMM_state,'Eta')
  # Eta = get_posterior_mean(Eta_samples)
  #
  # return(Eta[nas,1])
# }


  # results = data.frame(Method = 'MegaLMM',
  #                      pearson = cor(data$BLUP[nas],U[nas,1])/sqrt(h2_BLUE),
  #                      g_cor = calc_gcor(data,nas,U[nas,1],sKnn))
  # results$fold = foldid
  # write.csv(results,file = sprintf('%s/results_MegaLMM_fold_%d.csv',results_dir,foldid))


