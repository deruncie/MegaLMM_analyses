# This function takes observed trait values and genomic predictions and estimates the genetic correlation between them.
# This fulfills the parametric method for evaluation genomic prediction accuracy from Runcie and Cheng 2019

require(Matrix)
# require(rstan)

#' Title
#'
#' @param data data.frame with columns 'ID','obs','pred'
#' @param Knn kinship among the testing set, all data$ID must be in rownames(Knn)
#' @param sKnn optional; result of \code{svd(Knn)}
#' @param method analysis method to use
#' @param normalize Should \code{obs} and \code{pred} be standardized to unit variance?
#' @param control List of optional arguments to modify behavior of methods.
#'
#' @return
#' @export
#'
#' @examples
estimate_gcor = function(data,Knn,sKnn = NULL,method = c('MCMCglmm','sommer'),normalize = c(T,F),control = c()) {
  if(normalize) {
    data$obs = scale(data$obs)
    data$pred = scale(data$pred)
  }

  result = switch(method,
                  MCMCglmm = {
                    require(MCMCglmm)
                    if(is.null(sKnn)) sKnn = svd(Knn)

                    Dinv = as(diag(1/sKnn$d),'dgCMatrix')
                    rownames(Dinv) = colnames(Dinv) = rownames(Knn)
                    fixed = formula('cbind(ut_obs,ut_pred) ~ 0+ut1:trait')

                    n_test = nrow(data)
                    data$ut1 = t(sKnn$u) %*% matrix(1,nr = n_test,nc=1)
                    data$ut_obs = t(sKnn$u) %*% data$obs
                    data$ut_pred = t(sKnn$u) %*% data$pred

                    prior = list(R = list(V = diag(c(.5,.01),2),nu = 3),G = list(G1 = list(V = diag(c(.5,.5),2),nu = 3,alpha.mu = rep(0,2),alpha.V = diag(1,2))))#
                    nitt = 30000
                    thin = 50
                    burn = 5000
                    if('prior' %in% names(control)) prior = control$prior
                    if('nitt' %in% names(control)) nitt = control$nitt
                    if('thin' %in% names(control)) thin = control$thin
                    if('burn' %in% names(control)) burn = control$burn

                    m_gcor <- MCMCglmm(fixed,random = ~us(trait):ID,rcov = ~us(trait):units,data = data,prior = prior,family = rep('gaussian',2),pl=F,pr=T,
                                       ginverse = list(ID = Dinv),nitt = nitt,burn = burn,thin = thin,verbose = F)

                    h2_obs = apply(m_gcor$VCV,1,function(x) (x[1]/(x[1] + x[5])))
                    h2_pred = apply(m_gcor$VCV,1,function(x) (x[4]/(x[4] + x[8])))
                    g_cor = apply(m_gcor$VCV,1,function(x) x[2]/sqrt(x[1]*x[4]))
                    r_cor = apply(m_gcor$VCV,1,function(x) x[6]/sqrt(x[5]*x[8]))

                    return(c(g_cor = mean(g_cor * sqrt(h2_pred)), Rhat = rstan::Rhat(g_cor * sqrt(h2_pred))))

                  },
                  sommer = {
                    require(sommer)

                    data2 = rbind(data.frame(data,y = data$obs,Type = 'obs'),data.frame(data,y=data$pred,Type='pred'))
                    data2 = droplevels(data2)

                    g_cor3 = try({
                      ms = mmer(y ~ 0+Type,random = ~vs(us(Type),ID,Gu=Knn),rcov = ~vs(us(Type),units),data = data2,verbose = 0)
                      pin(ms,g_cor~V2/sqrt(V1*V3) * sqrt(V3/(V3+V6)))
                    },silent = T)
                    if(class(g_cor3) == 'try-error') {
                      ms = mmer(y ~ 0+Type,random = ~vs(us(Type),ID,Gu=Knn),rcov = ~vs(ds(Type),units),data = data2,verbose = 0)
                      g_cor3 = try(pin(ms,g_cor~V2/sqrt(V1*V3) * sqrt(V3/(V3+V5))),silent = T)
                      if(class(g_cor3) == 'try-error') g_cor3 = c(Estimate = NA,SE = NA)
                    }
                    return(unlist(g_cor3))

                  }
            )
  return(result)
}
