library(cowplot)
library(ggplot2)
library(dplyr)
library(emmeans)
library(car)
library(plyr)

results_dirs = c(MODEM = 'Results_MODEM', At_1001 = 'Results_At_1001')
# results_dirs = c(MODEM = 'Results_ARD_test_MODEM')
for(set in names(results_dirs)) {
  dir = results_dirs[set]

  results = do.call(bind_rows,lapply(list.files(pattern = 'Times',path=dir,full.names = T),function(x) read.csv(x)))
  results$Method = revalue(results$Method,c('BSFG'= 'MegaLMM'))
  results$Method = factor(results$Method,levels = c('MegaLMM','MCMCglmm','MTG2','phenix'))
  p_max = max(results$p)
  (p0 <- ggplot(results,aes(x=p,y=X2)) +
    scale_x_continuous(breaks = unique(results$p),trans = scales::log2_trans()) +
    geom_boxplot(aes(color=Method,group=interaction(Method,p)))
  )
  save_plot(sprintf('%s_ESS.pdf',set),p0)
  results$time2 = results$time
  results$time2[results$Method == 'MCMCglmm'] = with(subset(results,Method == 'MCMCglmm'),time*5/7 + time*2/7*1000/X2)
  results$time2[results$Method == 'MegaLMM'] = with(subset(results,Method == 'MegaLMM'),time*5/7 + time*2/7*1000/X2)

  results$log_time2 = log2(results$time2)
  results2 = aggregate(log_time2~Method+p,subset(results,Method != 'MegaLMM' & p %in% c(32,64)),FUN=mean)
  results2 = rbind(results2,data.frame(
    Method = c('MCMCglmm','MTG2'),p=512,
    log_time2 = c(predict(lm(log_time2~log2(p),subset(results2,Method == 'MCMCglmm')),newdata = list(p=(512))),
                  predict(lm(log_time2~log2(p),subset(results2,Method == 'MTG2')),newdata = list(p=(512))))
    # predict(lm(log_time2~log2(p),subset(results2,Method == 'phenix')),newdata = list(p=(512))))
  ))
  results3 = data.frame(p = 2^c(9:13,log2(p_max)));results3$time2 = results3$p/max(results3$p)*10*3600
  my_comma = function(x) {
    scales::comma(x,accuracy = 0.001,trim = F)
  }
  exponential_funs = expand.grid(p = 2^(10:13),slope = c(3.4,0.5,1,2,3))
  MegaLMM_1024 = exp(mean(log(subset(results,Method == 'MegaLMM' & p==exponential_funs$p[1])$time2)))
  # MegaLMM_1024 = mean(subset(results,Method == 'MegaLMM' & p==exponential_funs$p[1])$time2)
  exponential_funs$time2 = (exponential_funs$p/(2^10))^exponential_funs$slope*MegaLMM_1024
  (p1 <- ggplot(results,aes(x=p,y=time2/3600))  + scale_y_log10(breaks = 10^(-4:10),labels = my_comma) +
      scale_x_continuous(breaks = unique(results$p),trans = scales::log2_trans()) + ylab('Hours') + xlab('# traits') +
      # geom_line(aes(color=Method,group=Method),stat = 'summary',fun.y = 'mean',size=1)
      geom_smooth(aes(color=Method,group=Method),se=F) 
      +geom_line(data = results2,aes(x=p,y=(2^log_time2)/3600,color=Method,group=Method),linetype = 2,size=1)
      # +geom_line(data = results3,color = 'red',linetype=2)
    + geom_line(data = subset(exponential_funs,slope <= 3),aes(group = slope),color = 'grey70')
    + geom_text(data = subset(exponential_funs,p == 2^13 & slope <= 3),aes(x=1.3*p,y=time2/3600,label = slope))
    + geom_text(data = subset(exponential_funs,p == 2^13 & slope > 3),aes(x=1.3*p,y=time2/3600,label = 'Rate'))
      +geom_point(aes(color=Method,group=Method))
  )
  
  g_cors = do.call(bind_rows,lapply(list.files(pattern = 'Gcor',path=dir,full.names = T),function(x) cbind(x,read.csv(x))))
  g_cors$method = factor(g_cors$method,levels = c('BSFG','MCMCglmm','MTG2','phenix','rrBLUP'),labels = c('MegaLMM','MCMCglmm','MTG2','phenix','rrBLUP'))
  g_cors$g_cor3 = g_cors$MCMCglmm
  # g_cors$g_cor3 = g_cors$sommer
  # g_cors$g_cor3[g_cors$g_cor3>1] = NA

  #  position = position_dodge(width = .2)
  # (p2 <- ggplot(g_cors,aes(x=p,y=g_cor3)) +
  #     stat_summary(aes(group = method,color=method),fun.y=mean,geom='line',position = position) +
  #     stat_summary(aes(group = method,color=method),fun.ymin=function(x) mean(x)-sd(x)/sqrt(length(x))*2,fun.ymax=function(x) mean(x)+sd(x)/sqrt(length(x))*2,geom='errorbar',position = position) +
  #     scale_x_continuous(breaks = unique(g_cors$p),trans = scales::log2_trans())
  # )

  g_cors$trait = factor(g_cors$trait)
  g_cors$p2 = factor(g_cors$p)
  m = lm((g_cor3) ~ method:p2+p2:trait,subset(g_cors,!is.na(g_cor3) & !is.infinite(g_cor3) &  method != 'BSFG2' & !is.na(trait)))
  means = emmeans(m,specs = c('p2','method'))
  means_summary = as.data.frame(summary(means))
  means_summary$p = as.numeric(as.character(means_summary$p2))
  methods = unique(means_summary$method)
  
  position = position_dodge(width = .2)
  (p2 <- ggplot(subset(means_summary,method != 'rrBLUP' & !is.na(emmean)),aes(x=p,y=emmean)) + ylab('Estimated prediction accuracy') + xlab('# traits') +
      geom_errorbar(aes(color = method,ymin = lower.CL,ymax = upper.CL),position = position) +
      geom_point(aes(color = method),position = position) +
      geom_hline(yintercept=subset(means_summary,method == 'rrBLUP')$emmean[1]) +
      # geom_line(aes(color = method,group = method),position = position,size=1) +
      geom_smooth(aes(color=method,group=method),se=F) + 
      scale_x_continuous(breaks = unique(g_cors$p),trans = scales::log2_trans())
  )
  leg = get_legend(p1 + theme(legend.position = 'bottom'))
  nl = theme(legend.position = 'none')
  (p3 <- plot_grid(plot_grid(p2+nl,p1+nl,labels = c('A','B')),leg,nrow = 2,rel_heights = c(1,.2)))

  save_plot(filename = sprintf('%s_results.pdf',set),p3,base_asp = 2,base_height = 6)
}

