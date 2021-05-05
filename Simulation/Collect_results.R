library(cowplot)
library(ggplot2)
library(dplyr)
library(car)
library(plyr)
library(tidyr)
library(DescTools)


dir = 'Results_At_1001'
g_cors = do.call(bind_rows,lapply(list.files(pattern = 'Gcor_0',path=dir,full.names = T),function(x) cbind(x,read.csv(x))))
colnames(g_cors)[2] = 'Method'
write.csv(g_cors,file = 'G_cor_results.csv',row.names=F)

results = do.call(bind_rows,lapply(list.files(pattern = 'Times_At_1001_tmp_0',path=dir,full.names = T),function(x) tibble(File=x,read.csv(x))))
results$time2 = results$time
results$time2[results$Method == 'MCMCglmm'] = with(subset(results,Method == 'MCMCglmm'),time*5/7 + time*2/7*1000/X2)
results$time2[grepl('MegaLMM',results$Method)] = with(subset(results,grepl('MegaLMM',results$Method)),time*5/7 + time*2/7*1000/X2)
results = results[,c('File','Method','p','converged','time2')]
write.csv(results,file = 'Timing_results.csv',row.names=F)


dir = 'Results_G_accuracy_At_1001'

results = do.call(bind_rows,lapply(list.files(pattern = 'G_accuracy_2',path=dir,full.names = T),function(x) tibble(File=x,read.csv(x))))
results$Method = factor(results$Method,levels = unique(results$Method)[c(2,1,3,4)])
write.csv(results,file = 'G_accuracy_results.csv',row.names=F)

# Figure 1A: time comparisons among methods

results = read.csv('Timing_results.csv')
p_max = max(results$p)
# separate K results from base results
MegaLMM_results = droplevels(subset(results,grepl('MegaLMM',results$Method)))
MegaLMM_results$Method = sprintf('K_%0.2f',sapply(MegaLMM_results$Method,function(x) as.numeric(strsplit(as.character(x),'_')[[1]][2])))
results = droplevels(subset(results,Method %in% c('MTG2','MCMCglmm','phenix','MegaLMM_1.00_')))
results$Method = revalue(results$Method,c('MegaLMM_1.00_'= 'MegaLMM'))
results$Method = factor(results$Method,levels = c('MegaLMM','MCMCglmm','MTG2','phenix'))


results$log_time2 = log2(results$time2)
results2 = aggregate(log_time2~Method+p,subset(results,Method != 'MegaLMM' & p %in% c(32,64)),FUN=mean)
results2 = rbind(results2,data.frame(
  Method = c('MCMCglmm','MTG2'),p=512,
  log_time2 = c(predict(lm(log_time2~log2(p),subset(results2,Method == 'MCMCglmm')),newdata = list(p=(512))),
                predict(lm(log_time2~log2(p),subset(results2,Method == 'MTG2')),newdata = list(p=(512))))
))
results3 = data.frame(p = 2^c(9:13,log2(p_max)));results3$time2 = results3$p/max(results3$p)*10*3600
my_comma = function(x) {
  scales::comma(x,accuracy = 0.001,trim = F)
}
exponential_funs = expand.grid(p = 2^(10:13),slope = c(3.4,0.5,1,2,3))
MegaLMM_1024 = exp(mean(log(subset(results,Method == 'MegaLMM' & p==exponential_funs$p[1])$time2)))
exponential_funs$time2 = (exponential_funs$p/(2^10))^exponential_funs$slope*MegaLMM_1024
(p1 <- ggplot(results,aes(x=p,y=time2/3600))  + scale_y_log10(breaks = 10^(-4:10),labels = my_comma) +
    scale_x_continuous(breaks = unique(results$p),trans = scales::log2_trans()) + ylab('Hours') + xlab('# traits') +
    geom_smooth(aes(color=Method,group=Method),se=F) 
  +geom_line(data = results2,aes(x=p,y=(2^log_time2)/3600,color=Method,group=Method),linetype = 2,size=1)
  # +geom_line(data = results3,color = 'red',linetype=2)
  + geom_line(data = subset(exponential_funs,slope <= 3),aes(group = slope),color = 'grey70')
  + geom_text(data = subset(exponential_funs,p == 2^13 & slope <= 3),aes(x=1.3*p,y=time2/3600,label = slope))
  + geom_text(data = subset(exponential_funs,p == 2^13 & slope > 3),aes(x=1.3*p,y=time2/3600,label = 'Rate'))
  +geom_point(aes(color=Method,group=Method),position = position_dodge(width=.2)) + theme_cowplot() + background_grid(major = 'xy')
)

# Figure 1B accuracy comparisons among methods

g_cors = read.csv('G_cor_results.csv')

g_cors_MegaLMM = droplevels(subset(g_cors,grepl('MegaLMM',g_cors$Method)))
g_cors_MegaLMM$Method = sprintf('K_%0.2f',sapply(g_cors_MegaLMM$Method,function(x) as.numeric(strsplit(as.character(x),'_')[[1]][2])))
g_cors_rrBLUP = droplevels(subset(g_cors,Method %in% c('rrBLUP')))
g_cors = droplevels(subset(g_cors,Method %in% c('MTG2','MCMCglmm','phenix','MegaLMM_1.00_')))
g_cors$Method = revalue(g_cors$Method,c('MegaLMM_1.00_'= 'MegaLMM'))
g_cors$Method = factor(g_cors$Method,levels = c('MegaLMM','MCMCglmm','MTG2','phenix'))

# g_cors_wide = pivot_wider(subset(g_cors,Method == 'MegaLMM_1.00_'),id_cols = c('trait'),names_from='p',values_from = 'g_cor')

ks = data.frame(p = c(2^(2:13),20843))
ks$k = pmin(166,pmax(4,ks$p/2))
MegaLMM_means = aggregate(g_cor~p,subset(g_cors,Method == 'MegaLMM'),mean)
ks$height = cummax(MegaLMM_means$g_cor[match(ks$p,MegaLMM_means$p)] + .06)

dodge_width = .5
g_cors$p2 = factor(g_cors$p)
(p2 <- ggplot(g_cors,aes(x=log2(p),y=g_cor)) + ylab('Estimated prediction accuracy') + xlab('# traits') +
    geom_pointrange(fatten=2,stat="summary", fun.data="mean_se",aes(color = Method,group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
    geom_smooth(aes(color=Method,group=Method),se=F) +
    geom_pointrange(data = g_cors_rrBLUP, fatten=2,stat="summary", fun.data="mean_se",color = 'grey30',aes(group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
    geom_hline(yintercept=mean(g_cors_rrBLUP$g_cor)) +
    geom_text(data = ks,y=max(ks$height,na.rm=T),aes(label = k),color = scales::hue_pal()(4)[1]) + 
    scale_x_continuous(breaks = unique(log2(g_cors$p)),labels = unique(g_cors$p)) +
    theme_cowplot() + background_grid(major = 'xy')
)
# dev.off()

leg = get_legend(p2 + theme(legend.position = 'bottom',legend.justification = 'center'))
nl = theme(legend.position = 'none')
(p3 <- plot_grid(plot_grid(p2+nl,p1+nl,labels = c('A','B')),leg,nrow = 2,rel_heights = c(1,.2)))

set = 'At_1001'
save_plot(filename = sprintf('%s_results.pdf',set),p3,base_asp = 2,base_height = 7)


# Now with varying K
(p1b <- ggplot(MegaLMM_results,aes(x=p,y=time2/3600))  + scale_y_log10(breaks = 10^(-4:10),labels = my_comma) +
    scale_x_continuous(breaks = unique(MegaLMM_results$p),trans = scales::log2_trans()) + ylab('Hours') + xlab('# traits') +
    geom_smooth(aes(color=Method,group=Method),se=F) 
  +geom_point(aes(color=Method,group=Method),position = position_dodge(width=.4)) + theme_cowplot() + background_grid(major = 'xy')
)
# dev.off()


dodge_width = .5
g_cors$p2 = factor(g_cors$p)
(p2b <- ggplot(g_cors_MegaLMM,aes(x=log2(p),y=g_cor)) + ylab('Estimated prediction accuracy') + xlab('# traits') +
    geom_pointrange(fatten=2,stat="summary", fun.data="mean_se",aes(color = Method,group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
    geom_smooth(aes(color=Method,group=Method),se=F) +
    geom_hline(yintercept=mean(g_cors_rrBLUP$g_cor)) +
    scale_x_continuous(breaks = unique(log2(g_cors$p)),labels = unique(g_cors$p)) +
    theme_cowplot() + background_grid(major = 'xy')
)
# dev.off()

leg = get_legend(p2b + theme(legend.position = 'bottom',legend.justification = 'center'))
nl = theme(legend.position = 'none')
(p3b <- plot_grid(plot_grid(p2b+nl,p1b+nl,labels = c('A','B')),leg,nrow = 2,rel_heights = c(1,.2)))


set = 'At_1001'
save_plot(filename = sprintf('%s_changeK.pdf',set),p3b,base_asp = 2,base_height = 7)

full_g_cors = bind_rows(g_cors,subset(g_cors_MegaLMM,Method != 'K_1.00'))
full_g_cors$Method = factor(full_g_cors$Method,levels = unique(full_g_cors$Method)[c(2,3,4,1,5,6)])
(p4 <- ggplot(full_g_cors,aes(x=log2(p),y=FisherZ(g_cor))) +  scale_x_continuous(breaks = unique(log2(g_cors$p)),labels = unique(g_cors$p)) +
  geom_line(aes(group = trait,color=trait)) + facet_wrap(~Method) + 
  ylab('Estimated prediction accuracy (FisherZ-transformed)') + xlab('# traits') +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
)
save_plot(p4,file = 'At_1001_accuracy_by_gene_by_p.pdf',base_height = 5,base_asp = 1.7)

g_cors_wide = pivot_wider(bind_rows(g_cors,g_cors_rrBLUP),id_cols = c('p','trait'),names_from = 'Method',values_from = 'g_cor')
g_cors_tall = pivot_longer(g_cors_wide,cols = -c('p','trait','MegaLMM'))
g_cors_tall$Method = g_cors_tall$name
(p5 <- ggplot(subset(g_cors_tall,p < 1000),aes(y=MegaLMM,x=value)) + geom_point(aes(color = Method),size=1) + facet_wrap(~p) + geom_abline(slope=1,intercept=0) +
  scale_color_manual(values = c(scales::hue_pal()(4)[-1],'grey70'),drop=FALSE) + 
  ylab('MegaLMM') + xlab('Alternative method')# + theme(legend.position = c(.65,.1),legend.direction = 'horizontal')
  )

save_plot(plot_grid(p5+nl,get_legend(p5+theme(legend.position = 'bottom',legend.direction = 'horizontal')),ncol=1,rel_heights = c(1,.1)),
          file = 'At_1001_accuracy_by_gene_by_method.pdf',base_height = 5,base_asp = 1.7)
dev.off()


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
# leg = get_legend(p1 + theme(legend.position = 'bottom'))
# nl = theme(legend.position = 'none')
# (p3 <- plot_grid(plot_grid(p2+nl,p1+nl,labels = c('A','B')),leg,nrow = 2,rel_heights = c(1,.2)))
# 
# save_plot(filename = sprintf('%s_results.pdf',set),p3,base_asp = 2,base_height = 6)






results_dirs = c(MODEM = 'Results_MODEM', At_1001 = 'Results_At_1001')
# results_dirs = c(MODEM = 'Results_ARD_test_MODEM')
set = names(results_dirs)[2]
# for(set in names(results_dirs)) {
  dir = results_dirs[set]

  results = do.call(bind_rows,lapply(list.files(pattern = 'Times_At_1001_tmp_0',path=dir,full.names = T),function(x) tibble(File=x,read.csv(x))))
  results$Method = revalue(results$Method,c('BSFG'= 'MegaLMM'))
  # results$Method = factor(results$Method,levels = c('MegaLMM','MCMCglmm','MTG2','phenix'))
  p_max = max(results$p)
  (p0 <- ggplot(results,aes(x=p,y=X2)) +
    scale_x_continuous(breaks = unique(results$p),trans = scales::log2_trans()) +
    geom_boxplot(aes(color=Method,group=interaction(Method,p)))
  )
  save_plot(sprintf('%s_ESS.pdf',set),p0)
  results$time2 = results$time
  results$time2[results$Method == 'MCMCglmm'] = with(subset(results,Method == 'MCMCglmm'),time*5/7 + time*2/7*1000/X2)
  results$time2[grepl('MegaLMM',results$Method)] = with(subset(results,grepl('MegaLMM',results$Method)),time*5/7 + time*2/7*1000/X2)

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
    # + geom_line(data = subset(exponential_funs,slope <= 3),aes(group = slope),color = 'grey70')
    # + geom_text(data = subset(exponential_funs,p == 2^13 & slope <= 3),aes(x=1.3*p,y=time2/3600,label = slope))
    # + geom_text(data = subset(exponential_funs,p == 2^13 & slope > 3),aes(x=1.3*p,y=time2/3600,label = 'Rate'))
      +geom_point(aes(color=Method,group=Method)) + theme_cowplot() + background_grid(major = 'xy')
  )
  
  g_cors = do.call(bind_rows,lapply(list.files(pattern = 'Gcor_0',path=dir,full.names = T),function(x) cbind(x,read.csv(x))))
  # g_cors$method = factor(g_cors$method,levels = c('BSFG','MCMCglmm','MTG2','phenix','rrBLUP'),labels = c('MegaLMM','MCMCglmm','MTG2','phenix','rrBLUP'))
  g_cors$g_cor3 = g_cors$g_cor
  # g_cors$g_cor3 = g_cors$MCMCglmm
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
  # leg = get_legend(p1 + theme(legend.position = 'bottom'))
  # nl = theme(legend.position = 'none')
  # (p3 <- plot_grid(plot_grid(p2+nl,p1+nl,labels = c('A','B')),leg,nrow = 2,rel_heights = c(1,.2)))
  # 
  # save_plot(filename = sprintf('%s_results.pdf',set),p3,base_asp = 2,base_height = 6)
  
  library(multcomp)
  library(emmeans)
  library(DescTools)
  lapply(unique(g_cors$p),function(i) cld(emmeans(lm((g_cor3)~trait+method,subset(g_cors,p==i & grepl('MegaLMM',method))),specs = 'method')))
  
  g_cors_wide = pivot_wider(subset(g_cors,method != 'rrBLUP'),id_cols = c('p','trait'),names_from = 'method',values_from = 'g_cor3')
  g_cors_tall = pivot_longer(g_cors_wide,cols = -c('p','trait','MegaLMM_1.00_'))
  ggplot(subset(g_cors_tall,grepl('MegaLMM',name)),aes(x=p,y=MegaLMM_1.00_-value)) + geom_boxplot(aes(group = interaction(name,p),color=name),position = position_dodge(preserve = "single")) + scale_x_log10()
  ggplot(subset(g_cors_tall,grepl('MegaLMM',name)),aes(x=p,y=MegaLMM_1.00_-value)) + geom_jitter(aes(group = interaction(name,p),color=name),position = position_dodge()) + scale_x_log10()
  
  ggplot(g_cors_tall,aes(x=MegaLMM_1.00_,y=value)) + geom_point(aes(color = name)) + facet_wrap(~p) + geom_abline(slope=1,intercept=0)

  
  library(DescTools)
  library(numform)
  transform = function(x) FisherZ(x)
  
  g_cors_wide = pivot_wider(g_cors,id_cols = c('p','trait'),names_from = 'method',values_from = 'g_cor3')
  g_cors_tall = pivot_longer(g_cors_wide,cols = -c('p','trait','rrBLUP'))
  ggplot(g_cors_tall,aes(x=rrBLUP,y=value)) + geom_point(aes(color = name)) + facet_wrap(~p) + geom_abline(slope=1,intercept=0)
  
  
   g_cors_tall$gcor_effect = g_cors_tall$value - g_cors_tall$rrBLUP
  
  
  g_cors_tall$gcor_effect = g_cors_tall$value - g_cors_tall$rrBLUP
  
  
  g_cors_tall$trait = factor(g_cors_tall$trait)
  g_cors_tall$p2 = factor(g_cors_tall$p)
  g_cors_tall$method = g_cors_tall$name
  m = lm((gcor_effect) ~ method:p2+p2:trait,subset(g_cors_tall,!is.na(gcor_effect) & !is.infinite(gcor_effect) &  method != 'BSFG2' & !is.na(trait)))
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
      scale_x_continuous(breaks = unique(g_cors_tall$p),trans = scales::log2_trans())
  )
  
  g_cors$gcor_trans = FisherZ(g_cors$g_cor3)
  
  g_cors_wide = pivot_wider(g_cors,id_cols = c('p','trait'),names_from = 'method',values_from = 'gcor_trans')
  g_cors_tall = pivot_longer(g_cors_wide,cols = -c('p','trait','MegaLMM_1.00_'))
  g_cors_tall$gcor_effect = g_cors_tall$value - g_cors_tall$MegaLMM_1.00_
  
  
  g_cors_tall$trait = factor(g_cors_tall$trait)
  g_cors_tall$p2 = factor(g_cors_tall$p)
  g_cors_tall$method = g_cors_tall$name
  # m = lm((gcor_effect) ~ method:p2+p2:trait,subset(g_cors_tall,!is.na(gcor_effect) & !is.infinite(gcor_effect) &  method != 'BSFG2' & !is.na(trait)))
  # means = emmeans(m,specs = c('p2','method'))
  # means_summary = as.data.frame(summary(means))
  # means_summary$p = as.numeric(as.character(means_summary$p2))
  # methods = unique(means_summary$method)
  # 
  means_summary = aggregate(gcor_effect~method+p,data = subset(g_cors_tall,!is.na(gcor_effect)),function(x) {
    c(mean = mean(x),SE = sd(x)/sqrt(length(x)))
  })
  means_summary$mean = means_summary$gcor_effect[,'mean']
  means_summary$SE = means_summary$gcor_effect[,'SE']
  means_summary$lower.CL = means_summary$mean - 2*means_summary$SE
  means_summary$upper.CL = means_summary$mean + 2*means_summary$SE
  means_summary$p = as.numeric(as.character(means_summary$p))
  
  position = position_dodge(width = .2)
  (p2 <- ggplot(means_summary,aes(x=p,y=mean)) + ylab('Estimated prediction accuracy') + xlab('# traits') +
      geom_errorbar(aes(color = method,ymin = lower.CL,ymax = upper.CL),position = position) +
      geom_point(aes(color = method),position = position) +
      # geom_hline(yintercept=subset(means_summary,method == 'rrBLUP')$emmean[1]) +
      # geom_line(aes(color = method,group = method),position = position,size=1) +
      geom_smooth(aes(color=method,group=method),se=F) +
      scale_x_continuous(breaks = unique(g_cors_tall$p),trans = scales::log2_trans())
  )
  g_cors_bygene = pivot_wider(g_cors,id_cols = c('trait','method'),names_from = 'p',values_from = 'gcor_trans')
  g_cors_bygene_tall = pivot_longer(g_cors_bygene,cols = -c(1:3))
  
  g_cors_bygene_tall$rel_value = g_cors_bygene_tall$value# - g_cors_bygene_tall$`4`
  g_cors_bygene_tall$p = as.numeric(as.character(g_cors_bygene_tall$name))
  ggplot(g_cors_bygene_tall,aes(x = p,y=rel_value)) + scale_x_log10() + geom_line(aes(group = trait)) + facet_wrap(~method)
  
  
  g_cors_wide = pivot_wider(g_cors,id_cols = c('p','trait'),names_from = 'method',values_from = 'gcor_trans')
  g_cors_tall = pivot_longer(g_cors_wide,cols = -c('p','trait','rrBLUP'))
  g_cors_tall$gcor_norm = g_cors_tall$value - g_cors_tall$rrBLUP
  ggplot(g_cors_tall,aes(x = p,y=gcor_norm)) + scale_x_log10() + geom_line(aes(group = trait)) + facet_wrap(~name)
  
  ggplot(subset(g_cors,method != 'rrBLUP' & trait %in% unique(g_cors$trait)[1:8]),aes(x=p,y=g_cor3)) + geom_line(aes(group = method,color = method)) + scale_x_log10() + facet_wrap(~trait)
  
   # }


(p1 <- ggplot(subset(results,Method %in% c('MCMCglmm','MTG2','MegaLMM')),aes(x=p,y=time2/3600))  + 
    scale_y_log10(breaks = 10^(-4:4),labels = my_comma,limits = c(10^-4,10^4)) + 
    coord_cartesian(ylim=c(1e-4,1e2)) + 
    scale_x_continuous(breaks = unique(results$p),trans = scales::log2_trans(),limits = range(results$p)) + ylab('Hours') + xlab('# traits') +
    # geom_line(aes(color=Method,group=Method),stat = 'summary',fun.y = 'mean',size=1)
    geom_smooth(aes(color=Method,group=Method),se=F) + 
    scale_color_discrete(drop=F) 
  +geom_line(data = subset(results2,Method != 'phenix'),aes(x=p,y=(2^log_time2)/3600,color=Method,group=Method),linetype = 2,size=1)
  # +geom_line(data = results3,color = 'red',linetype=2)
  # + geom_line(data = subset(exponential_funs,slope <= 3),aes(group = slope),color = 'grey70')
  # + geom_text(data = subset(exponential_funs,p == 2^13 & slope <= 3),aes(x=1.3*p,y=time2/3600,label = slope))
  # + geom_text(data = subset(exponential_funs,p == 2^13 & slope > 3),aes(x=1.3*p,y=time2/3600,label = 'Rate'))
  +geom_point(aes(color=Method,group=Method)) + theme_cowplot() + background_grid(major = 'xy')
)

#,'MCMCglmm','phenix','MegaLMM'
rrBLUP = subset(means_summary,method == 'rrBLUP')
means_summary = droplevels(subset(means_summary,method != 'rrBLUP'))
(p2 <- ggplot(subset(means_summary,method %in% c('MTG2','MCMCglmm') & !is.na(emmean)),aes(x=p,y=emmean)) + ylab('Estimated prediction accuracy') + xlab('# traits') +
    geom_errorbar(aes(color = method,ymin = lower.CL,ymax = upper.CL),position = position) +
    geom_point(aes(color = method),position = position) +
    geom_hline(yintercept=subset(rrBLUP,method == 'rrBLUP')$emmean[1]) +
    ylim(c(0.2,.8)) + 
    scale_x_continuous(breaks = unique(means_summary$p),trans = scales::log2_trans()) + coord_cartesian(xlim=c(3,66)) + 
    geom_line(aes(color = method,group = method),position = position,size=1) +
    # geom_smooth(aes(color=method,group=method),se=F) + 
    scale_color_discrete(drop=F) + 
    theme_cowplot()+ background_grid(major = 'xy') + theme(legend.position = 'none')
)

library(tidyverse)
g_cors_wide = pivot_wider(g_cors,names_from = 'method',id_cols = c('p','trait'),values_from = 'g_cor')


results_wide = pivot_wider(results,id_cols = c('p','File'),names_from = 'Method',values_from = 'log_time2')
results_comp = pivot_longer(results_wide,cols = -c('p','File','MegaLMM_1.00_'))
ggplot(results_comp,aes(x=MegaLMM_1.00_,y = value)) +geom_point(aes(color = p)) + geom_abline(slope=1,intercept=0) + facet_wrap(~name,scales= 'free')


g_cors = do.call(bind_rows,lapply(list.files(pattern = 'Gcor_prior',path=dir,full.names = T),function(x) cbind(x,read.csv(x))))

g_cors = pivot_wider(g_cors,id_cols = c('p','k','trait'),names_from = c('prop_0','delta_scale'),values_from = 'g_cor')
g_cors = data.frame(g_cors)
g33 = subset(g_cors,p==33)
g33[,4:7] = g33[,4:7]/g33[,5]
g8193 = subset(g_cors,p==8193)
g8193[,4:7] = g8193[,4:7]/g8193[,5]
