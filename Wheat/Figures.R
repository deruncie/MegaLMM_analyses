library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpattern)

# Main figure
# A: Performance of different methods in 1 trial
# B: Gcor and Pcor
results = read.csv(file = 'Results/Results_1415_OF.csv')
results = subset(results,Method != 'MegaLMM\nHorseshoe')
results$Method = factor(results$Method,levels = unique(results$Method)[c(3,8,7,4,5,1,2,6)])
results$HTP = results$Method %in% levels(results$Method)[-c(1:3)]
results$CV1 = results$Method %in% levels(results$Method)[8]
results$MegaLMM = results$Method %in% levels(results$Method)[-c(6:8)]

# library(DescTools)
(p1 <- ggplot(results,aes(x=Method,y=g_cor)) + 
    # stat_summary(geom = 'bar',fun = 'mean',aes(group = interaction(HTP,MegaLMM,Method),fill = !MegaLMM,color = CV1)) +
    # geom_bar(stat = 'summary',fun.y = 'mean',aes(group = interaction(HTP,MegaLMM,Method),fill = !MegaLMM)) +
    geom_bar_pattern(stat = 'summary',fun = 'mean',aes(group = interaction(HTP,MegaLMM,Method),fill = MegaLMM,pattern = CV1),
                     pattern_color = "grey40",
                     pattern_angle = 35,
                     pattern_density = 0.01,
                     pattern_spacing = 0.02) +
    geom_errorbar(stat="summary", fun.data="mean_se",width=.3) +
    facet_grid(~HTP,scales = 'free_x',space = 'free_x',labeller = labeller(HTP = c(`TRUE` = 'Genomic data \nplus hyperspectral data',`FALSE` = 'Genomic data \nonly'))) +
    guides(fill = F,pattern=F) + #ylim(c(0,1)) +
    # scale_color_manual(values = c(NA,'black')) + 
    scale_pattern_manual(values = c('none','stripe')) + 
    scale_fill_manual(values = scales::brewer_pal(palette = 'Set2')(4)[c(2,1)],drop=FALSE) + 
    # xlab('Method') +
    xlab('') +
    ylab('Estimated prediction accuracy')+ theme_bw() + 
    theme(axis.text.x = element_text(angle=45,vjust = .5),panel.grid.major = element_line())) 

G_cor_results = read.csv('G_cor_results.csv')
(p2 <- ggplot(G_cor_results,aes(x=wavelength)) +
    facet_wrap(~date) +
    geom_hline(yintercept = 0,size=.25) +
    geom_ribbon(aes(ymin = G_cor_low,ymax = G_cor_high),alpha = 0.3) +
    geom_line(aes(y = G_cor_mean,col = 'Genetic')) +
    geom_line(aes(y = P_cor,col = 'Phenotypic')) +
    scale_color_manual(values=c('Genetic' = 'red','Phenotypic' = 'black'),name = 'Correlation') +
    xlab('Wavelength') + ylab('Correlation') +
    theme_bw() + theme(legend.position = c(.8,.15)))

(p <- plot_grid(p1,p2,labels = c('A','B')))
save_plot(p,file = 'Figures/Krause_results_panels.pdf',base_asp = 2,base_height = 7)


# Supplements: all trials
results = read.csv('Results/All_trials_results.csv')
results = subset(results,fold < 6)
results$Method = factor(results$Method,levels = unique(results$Method)[c(1,8,7,2,3,4,6,5)])
results$Method2 = factor(results$Method,levels = levels(results$Method),labels = sub('\n',' ',levels(results$Method)))
results$HTP = results$Method %in% levels(results$Method)[-c(1:3)]
results$CV1 = results$Method %in% levels(results$Method)[8]
results$MegaLMM = results$Method %in% levels(results$Method)[-c(6:8)]

(p3 <- ggplot(results,aes(x=Method2,y=g_cor)) + 
  # geom_boxplot(aes(group = interaction(Method2,Trial),color = Method2)) + 
  geom_bar_pattern(stat = 'summary',fun = 'mean',aes(group = interaction(HTP,MegaLMM,Method),fill = MegaLMM,pattern = CV1),
  # geom_boxplot_pattern(aes(group = interaction(HTP,MegaLMM,Method),fill = !MegaLMM,pattern = CV1),
                   pattern_color = "grey40",
                   pattern_angle = 35,
                   pattern_density = 0.01,
                   pattern_spacing = 0.05) +
  geom_errorbar(stat="summary", fun.data="mean_se",width=.3) +
  theme(axis.text.x = element_text(angle=90,hjust = 1),legend.position = 'none') +
  facet_grid(trt~year,scales = 'free') + xlab('') + ylab('Estimated prediction accuracy') + labs(color = 'Method')+
  guides(fill = F,pattern=F) + #ylim(c(0,1)) +
  scale_pattern_manual(values = c('none','stripe')) + 
  scale_fill_manual(values = scales::brewer_pal(palette = 'Set2')(4)[c(2,1)],drop=FALSE) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90,hjust=1)))
save_plot(p3,file = 'Figures/Krause_allTrials.pdf',base_height=8,base_asp=1)
