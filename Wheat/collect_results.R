library(data.table)
library(foreach)
library(emmeans)
library(tidyr)
library(ggplot2)
library(plyr)
library(cowplot)

# load results for 1415_OF
results_folder = 'Results_1415_OF_Bgcor'

results = foreach(file = list.files(path=results_folder),.combine = 'rbind') %do% {
  a = fread(sprintf('%s/%s',results_folder,file))
  if(grepl('ARD',file)) a$Method = paste(a$Method,'ARD',sep='_')
  a
}

results = subset(results,Method != 'GBLUP_rrBLUP')
results$Method = revalue(results$Method,c('BSFG' = 'MegaLMM\nGBLUP',
                                          'BSFG_RKHS' = 'MegaLMM\nRKHS',
                                          'BSFG_X' = 'MegaLMM\nHorseshoe',
                                          'GBLUP_H' = 'HBLUP',
                                          'GBLUP_KH' = 'GBLUP+H',
                                          'MegaLMM_K_0' = 'MegaLMM\nGBLUP CV1',
                                          'BL' = 'Bayesian\nLasso')) #,'GBLUP_Hbest' = 'GBLUP_H6'))
results = subset(results,Method %in% c('GBLUP_KHbest', 'GBLUP_Hbest') == F)
write.csv(results[,-1],file = 'Results_1415_OF.csv',row.names = F)

# load correlations from full-data run
# files are in /home/deruncie/projects/BSFG/Krause
source('data_prep_Krause.R')
G_cor_samples = readRDS('BSFG_Krause_K_full_1/G_cor_samples.rds')
G_cor_mean = colMeans(G_cor_samples)
G_cor_HPD = get_posterior_HPDinterval(G_cor_samples)
P_cor = cor(data$BLUE,HTP_wide[,-1])[1,]
results = data.frame(trait = names(P_cor),G_cor_low = G_cor_HPD[1,-1],G_cor_high = G_cor_HPD[2,-1],G_cor_mean = G_cor_mean[-1],P_cor = P_cor)
results = separate(results,'trait',c('wavelength','date'),sep='::')
results$wavelength = as.numeric(substr(results$wavelength,12,14))
results$date = factor(results$date,levels = unique(results$date),labels = format(as.Date(unique(results$date),format = '%y%m%d'),'%d-%b'))
write.csv(results,file = '/group/runciegrp2/Projects/MegaLMM_revisions/Wheat/G_cor_results.csv',row.names=F)



# load results from all trials

BLUEs = fread('Data/Krause_et_al_2018_Yield_BLUEs.csv',data.table=F)
BLUEs$Trial = paste(BLUEs$`Breeding Cycle`,BLUEs$Managed_Treatment,sep='\n')
trial_names = unique(BLUEs$Trial)
names(trial_names) = sprintf('./Trial_%02d',1:length(trial_names))
years = sapply(trial_names,function(x) strsplit(x,'\n')[[1]][1])
names(years) = sprintf('./Trial_%02d',1:length(years))
trts = sapply(trial_names,function(x) strsplit(x,'\n')[[1]][2])
names(trts) = sprintf('./Trial_%02d',1:length(trts))

results = c()
trials = list.dirs(recursive = F)
trials = trials[grep('Trial',trials)]
for(trial in trials) {
  files = list.files(path = sprintf('%s/Results',trial),pattern = '.csv',full.names = T)
  for(file in files) {
    name = strsplit(file,'/',fixed=T)[[1]];name = name[length(name)];name = sub('results_','',name);name = sub('_fold_1.csv','',name)
    results = rbind(results,data.frame(Trial = trial,Name = name, read.csv(file,row.names = 1)))
  }
}
results$year = years[results$Trial]
results$trt = trts[results$Trial]
results$Trial = trial_names[results$Trial]
results = droplevels(subset(results,Method != 'GBLUP_rrBLUP'))
results$Method = revalue(results$Method,c('MegaLMM_K' = 'MegaLMM\nGBLUP',
                                          'MegaLMM_K_0' = 'MegaLMM\nGBLUP CV1',
                                          'MegaLMM_RKHS' = 'MegaLMM\nRKHS',
                                          'GBLUP_H' = 'HBLUP',
                                          'GBLUP_KH' = 'GBLUP+H',
                                          'BL' = 'Bayesian\nLasso')) #,
results$Method = factor(results$Method,levels = levels(results$Method)[c(1,7,8,5,2,3,4,6)])

results = write.csv(results,file = 'All_trials_results.csv',quote=T,row.names=F)




# g_wide = pivot_wider(results,id_cols = 'fold',names_from = 'Method',values_from = 'g_cor')
# 
# g_tall = pivot_longer(g_wide[,-1],cols = everything())
# g_tall = subset(g_tall,!grepl('best',name))
# g_tall$name = revalue(g_tall$name,c('BSFG' = 'MegaLMM\nGBLUP','BSFG_X' = 'MegaLMM\nHorseshoe','BSFG_RKHS' = 'MegaLMM\nRKHS','GBLUP_H' = 'HBLUP','GBLUP_KH' = 'GBLUP+H')) #,'GBLUP_Hbest' = 'GBLUP_H6'
# g_tall$name = factor(g_tall$name,levels = unique(g_tall$name)[c(4,8,7,5,6,1,2,3)])
# g_tall$HTP = g_tall$name %in% c('MegaLMM\nGBLUP','MegaLMM\nHorseshoe','MegaLMM\nRKHS','HBLUP','GBLUP+H')
# g_tall$MegaLMM = g_tall$name %in% c('MegaLMM\nGBLUP','MegaLMM\nHorseshoe','MegaLMM\nRKHS')
# 
# (p1 <- ggplot(na.omit(g_tall),aes(x=name,y=value)) + geom_bar(stat = 'summary',fun.y = 'mean',aes(group = interaction(HTP,MegaLMM,name),fill = !MegaLMM)) +
#     geom_errorbar(stat="summary", fun.data="mean_se",width=.3) +
#     facet_grid(~HTP,scales = 'free_x',space = 'free_x',labeller = labeller(HTP = c(`TRUE` = 'Genomic data \nplus hyperspectral data',`FALSE` = 'Genomic data \nonly'))) +
#     guides(fill = F) + #ylim(c(0,1)) +
#     # xlab('Method') +
#     xlab('') +
#     ylab('Estimated prediction accuracy')) + theme(axis.text.x = element_text(angle=45,vjust = .5))
# 
# 
# 
# 
# dirs = list.dirs()
# dirs = dirs[grep('Posterior',dirs)]
# vars = c()
# for(dir in dirs) {
#   print(match(dir,dirs)/length(dirs))
#   files = list.files(path = dir,pattern = 'Lambda_',full.names=T)
#   Lambda = MegaLMM::get_posterior_mean(readRDS(tail(files,n=1)))
#   v = rowSums(Lambda^2) / sum(Lambda^2)
#   vars = rbind(vars,data.frame(dir=dir,K = 1:nrow(Lambda),perc_var = v,cumperc = cumsum(v)))
# }
# vars = tidyr::separate(vars,dir,c('Base','Trial','Method','Posterior'),sep='/',remove = F)
# vars$Method = sub('_[1-5]','',vars$Method)
# ggplot(vars,aes(x=K,y=perc_var)) + geom_line(aes(group = dir,color = Method)) + facet_wrap(~Trial) + scale_y_log10()
# dev.off()
# 
# 
# # results_gcor = as.data.frame(pivot_wider(subset(results,Method != 'MegaLMM'),id_cols = c('Trial','fold'),names_from = 'Method',values_from = 'g_cor'))
# library(ggplot2)
# # ggplot(results,aes(x=g_cor,y=pearson)) + geom_point() + geom_abline(slope=1,intercept = 0) + facet_wrap(~Method)
# # results_gcor_tall = as.data.frame(pivot_longer(results_gcor,cols = -c('Trial','fold','MegaLMM\nGBLUP')))
# # ggplot(results_gcor_tall,aes(x=`MegaLMM\nGBLUP`,y=value)) + geom_point() + geom_abline(slope=1,intercept = 0) + facet_wrap(~name)
# 
# # ggplot(results,aes(x=Method,y=g_cor)) + geom_boxplot(aes(group = interaction(Method,Trial),color = Method)) + facet_wrap(~Trial,scales = 'free')
# results$Method2 = factor(results$Method,levels = levels(results$Method),labels = sub('\n',' ',levels(results$Method)))
# ggplot(results,aes(x=Method2,y=g_cor)) + geom_boxplot(aes(group = interaction(Method2,Trial),color = Method2)) + 
#   theme(axis.text.x = element_text(angle=90,hjust = 1),legend.position = 'none') +
#   facet_grid(trt~year,scales = 'free') + xlab('') + ylab('Estimated prediction accuracy') + labs(color = 'Method')
# 
# 
# results_ave = aggregate(g_cor~Trial+Method,data = results,FUN = mean)
# ave_gcor_wide = as.data.frame(pivot_wider(results_ave,id_cols = c('Trial'),names_from = 'Method',values_from = 'g_cor'))
# ave_gcor_tall_comp = as.data.frame(pivot_longer(ave_gcor_wide,cols = -c('Trial','MegaLMM\nGBLUP')))
# ave_gcor_tall_comp$name = factor(ave_gcor_tall_comp$name,levels = unique(ave_gcor_tall_comp$name)[c(1,2,3,4,5,6,7)])
# ggplot(ave_gcor_tall_comp,aes(x=`MegaLMM\nGBLUP`,y=value)) + geom_point() + geom_abline(slope=1,intercept = 0) + facet_wrap(~name)
# 
# 
# results_folder = 'Results_1415_OF_Bgcor'
# 
# results = foreach(file = list.files(path=results_folder),.combine = 'rbind') %do% {
#   a = fread(sprintf('%s/%s',results_folder,file))
#   if(grepl('ARD',file)) a$Method = paste(a$Method,'ARD',sep='_')
#   a
# }
# 
# results = subset(results,Method != 'GBLUP_rrBLUP')
# 
# results$Method = sub('_ARD','',results$Method)
# g_wide = pivot_wider(results,id_cols = 'fold',names_from = 'Method',values_from = 'g_cor')
# 
# g_tall = pivot_longer(g_wide[,-1],cols = everything())
# g_tall = subset(g_tall,!grepl('best',name))
# g_tall$name = revalue(g_tall$name,c('BSFG' = 'MegaLMM\nGBLUP','BSFG_X' = 'MegaLMM\nHorseshoe','BSFG_RKHS' = 'MegaLMM\nRKHS','GBLUP_H' = 'HBLUP','GBLUP_KH' = 'GBLUP+H')) #,'GBLUP_Hbest' = 'GBLUP_H6'
# g_tall$name = factor(g_tall$name,levels = unique(g_tall$name)[c(4,8,7,5,6,1,2,3)])
# g_tall$HTP = g_tall$name %in% c('MegaLMM\nGBLUP','MegaLMM\nHorseshoe','MegaLMM\nRKHS','HBLUP','GBLUP+H')
# g_tall$MegaLMM = g_tall$name %in% c('MegaLMM\nGBLUP','MegaLMM\nHorseshoe','MegaLMM\nRKHS')
# 
# (p1 <- ggplot(na.omit(g_tall),aes(x=name,y=value)) + geom_bar(stat = 'summary',fun.y = 'mean',aes(group = interaction(HTP,MegaLMM,name),fill = !MegaLMM)) +
#   geom_errorbar(stat="summary", fun.data="mean_se",width=.3) +
#   facet_grid(~HTP,scales = 'free_x',space = 'free_x',labeller = labeller(HTP = c(`TRUE` = 'Genomic data \nplus hyperspectral data',`FALSE` = 'Genomic data \nonly'))) +
#   guides(fill = F) + #ylim(c(0,1)) +
#   # xlab('Method') +
#   xlab('') +
#   ylab('Estimated prediction accuracy')) + theme(axis.text.x = element_text(angle=45,vjust = .5))
# 
# # save_plot(p1,file = 'Krause_performance.pdf')
# aggregate(value~name,g_tall,FUN=mean)
# 
# source('data_prep_Krause.R')
# G_cor_samples = readRDS('BSFG_Krause_K_full_1/G_cor_samples.rds')
# G_cor_mean = colMeans(G_cor_samples)
# G_cor_HPD = get_posterior_HPDinterval(G_cor_samples)
# P_cor = cor(data$BLUE,HTP_wide[,-1])[1,]
# results = data.frame(trait = names(P_cor),G_cor_low = G_cor_HPD[1,-1],G_cor_high = G_cor_HPD[2,-1],G_cor_mean = G_cor_mean[-1],P_cor = P_cor)
# results = separate(results,'trait',c('wavelength','date'),sep='::')
# results$wavelength = as.numeric(substr(results$wavelength,12,14))
# results$date = factor(results$date,levels = unique(results$date),labels = format(as.Date(unique(results$date),format = '%y%m%d'),'%d-%b'))
# (p2 <- ggplot(results,aes(x=wavelength)) +
#     facet_wrap(~date) +
#     geom_hline(yintercept = 0,size=.25) +
#     geom_ribbon(aes(ymin = G_cor_low,ymax = G_cor_high),alpha = 0.3) +
#     geom_line(aes(y = G_cor_mean,col = 'Genetic')) +
#     geom_line(aes(y = P_cor,col = 'Phenotypic')) +
#     scale_color_manual(values=c('Genetic' = 'red','Phenotypic' = 'black'),name = 'Correlation') +
#     xlab('Wavelength') + ylab('Correlation') +
#     theme(legend.position = c(.8,.15)))
# 
# (p <- plot_grid(p1,p2,labels = c('A','B')))
# save_plot(p,file = 'Figures/Krause_results_panels.pdf',base_asp = 2,base_height = 7)
