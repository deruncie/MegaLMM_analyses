
library(cowplot)
library(ggplot2)
library(plyr)
library(dplyr)
library(DescTools)
library(tidyr)



# Figure 1A: time comparisons among methods

results = read.csv('Results/Timing_results.csv')
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

g_cors = read.csv('Results/G_cor_results.csv')

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


results = read.csv('Results/G_accuracy_results.csv')

dodge_width = .5
(p1 <- ggplot(results,aes(x=log2(p),y=G_RMSE_G)) + ylab('RMSE of genetic covariances') + xlab('# traits') + #ggtitle('G') + 
    # geom_point(aes(color = Method,group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
    geom_pointrange(fatten=2,stat="summary", fun.data="mean_se",aes(color = Method,group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
    geom_smooth(aes(color=Method,group=Method),se=F) +
    scale_x_continuous(breaks = unique(log2(results$p)),labels = unique(results$p)) +
    theme_cowplot() + background_grid(major = 'xy')
)

(p2 <- ggplot(results,aes(x=log2(p),y=R_RMSE_G)) + ylab('RMSE of residual covariances') + xlab('# traits')  + #ggtitle('R') + 
    geom_pointrange(fatten=2,stat="summary", fun.data="mean_se",aes(color = Method,group=Method),position = position_dodge2(width=dodge_width,preserve = 'single')) +
    geom_smooth(aes(color=Method,group=Method),se=F) +
    scale_x_continuous(breaks = unique(log2(results$p)),labels = unique(results$p)) +
    theme_cowplot() + background_grid(major = 'xy')
)
leg = get_legend(p1 + theme(legend.position = 'bottom',legend.justification = 'center'))
nl = theme(legend.position = 'none')
save_plot(plot_grid(plot_grid(p1+nl,p2+nl,labels = c('A','B')),leg,nrow = 2,rel_heights = c(1,.2)),file = 'Figures/G_accuracy.pdf',base_height = 5)
