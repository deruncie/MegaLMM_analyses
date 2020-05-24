library(data.table)
library(foreach)
library(emmeans)
library(tidyr)
library(ggplot2)
library(plyr)
library(cowplot)

results_folder = 'Results_1415_OF_Bgcor'

results = foreach(file = list.files(path=results_folder),.combine = 'rbind') %do% {
  a = fread(sprintf('%s/%s',results_folder,file))
  if(grepl('ARD',file)) a$Method = paste(a$Method,'ARD',sep='_')
  a
}

results = subset(results,Method != 'GBLUP_rrBLUP')

results$Method = sub('_ARD','',results$Method)
g_wide = pivot_wider(results,id_cols = 'fold',names_from = 'Method',values_from = 'g_cor')

g_tall = pivot_longer(g_wide[,-1],cols = everything())
g_tall = subset(g_tall,!grepl('best',name))
g_tall$name = revalue(g_tall$name,c('BSFG' = 'MegaLMM\nGBLUP','BSFG_X' = 'MegaLMM\nHorseshoe','BSFG_RKHS' = 'MegaLMM\nRKHS','GBLUP_H' = 'HBLUP','GBLUP_KH' = 'GBLUP+H')) #,'GBLUP_Hbest' = 'GBLUP_H6'
g_tall$name = factor(g_tall$name,levels = unique(g_tall$name)[c(4,8,7,5,6,1,2,3)])
g_tall$HTP = g_tall$name %in% c('MegaLMM\nGBLUP','MegaLMM\nHorseshoe','MegaLMM\nRKHS','HBLUP','GBLUP+H')
g_tall$MegaLMM = g_tall$name %in% c('MegaLMM\nGBLUP','MegaLMM\nHorseshoe','MegaLMM\nRKHS')

(p1 <- ggplot(na.omit(g_tall),aes(x=name,y=value)) + geom_bar(stat = 'summary',fun.y = 'mean',aes(group = interaction(HTP,MegaLMM,name),fill = !MegaLMM)) +
  geom_errorbar(stat="summary", fun.data="mean_se",width=.3) +
  facet_grid(~HTP,scales = 'free_x',space = 'free_x',labeller = labeller(HTP = c(`TRUE` = 'Genomic data \nplus hyperspectral data',`FALSE` = 'Genomic data \nonly'))) +
  guides(fill = F) + #ylim(c(0,1)) +
  # xlab('Method') +
  xlab('') +
  ylab('Estimated prediction accuracy')) + theme(axis.text.x = element_text(angle=45,vjust = .5))

# save_plot(p1,file = 'Krause_performance.pdf')
aggregate(value~name,g_tall,FUN=mean)

source('data_prep_Krause.R')
G_cor_samples = readRDS('BSFG_Krause_K_full_1/G_cor_samples.rds')
G_cor_mean = colMeans(G_cor_samples)
G_cor_HPD = get_posterior_HPDinterval(G_cor_samples)
P_cor = cor(data$BLUE,HTP_wide[,-1])[1,]
results = data.frame(trait = names(P_cor),G_cor_low = G_cor_HPD[1,-1],G_cor_high = G_cor_HPD[2,-1],G_cor_mean = G_cor_mean[-1],P_cor = P_cor)
results = separate(results,'trait',c('wavelength','date'),sep='::')
results$wavelength = as.numeric(substr(results$wavelength,12,14))
results$date = factor(results$date,levels = unique(results$date),labels = format(as.Date(unique(results$date),format = '%y%m%d'),'%d-%b'))
(p2 <- ggplot(results,aes(x=wavelength)) +
    facet_wrap(~date) +
    geom_hline(yintercept = 0,size=.25) +
    geom_ribbon(aes(ymin = G_cor_low,ymax = G_cor_high),alpha = 0.3) +
    geom_line(aes(y = G_cor_mean,col = 'Genetic')) +
    geom_line(aes(y = P_cor,col = 'Phenotypic')) +
    scale_color_manual(values=c('Genetic' = 'red','Phenotypic' = 'black'),name = 'Correlation') +
    xlab('Wavelength') + ylab('Correlation') +
    theme(legend.position = c(.8,.15)))

(p <- plot_grid(p1,p2,labels = c('A','B')))
save_plot(p,file = 'Figures/Krause_results_panels.pdf',base_asp = 2,base_height = 7)
