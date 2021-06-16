library(foreach)
library(doParallel)
library(ggplot2)
library(mapdata)
library(data.table)
library(usmap)
library(cowplot)
library(ggthemes)
library(dplyr)
library(tidyr)
source('my_pylo2map.R')
library(lme4)
library(lme4qtl)
library(geosphere)
registerDoParallel(RcppParallel::defaultNumThreads()-1)

files = list.files(path='BLUP_matrices/',pattern = 'csv')
traits = sub('_BLUPs.csv','',files,fixed=T)
traits = gsub(' ','_',traits)
traits = traits[c(4,1,2,3)]
trait_names = traits
traits = c('DTS','ASI','Grain Yield','Plant Height')
names(trait_names) = traits

a=foreach(trait. = trait_names) %do% {
  print(trait.)
  print(dim(fread(sprintf('BLUP_matrices/%s_BLUPs.csv',gsub('_',' ',trait.)))))
}

trait. = traits[1]
results = foreach(trait. = traits,.combine = 'bind_rows') %dopar% {
  trait_name = trait_names[trait.]
  foreach(i=1:20,.combine = 'bind_rows',.errorhandling = 'pass') %do% {
    try({
      res_BSFG = readRDS(sprintf('G2F_byTrait_results/%s_%d.rds',trait_name,i))
      res_alt = readRDS(sprintf('G2F_byTrait_results/alternatives_%s_%d.rds',trait_name,i))
      res_mean = readRDS(sprintf('G2F_byTrait_results/mean_%s_%d.rds',trait_name,i))

      Y = res_BSFG$Y
      mask = res_BSFG$mask
      Y[!mask] = NA
      BSFG_Eta = sapply(1:ncol(Y),function(i) cor(Y[,i],res_BSFG$Eta_m[,i],use='p'))
      BSFG_U = sapply(1:ncol(Y),function(i) cor(Y[,i],res_BSFG$U[,i],use='p'))
      rrBLUP = sapply(1:ncol(Y),function(i) cor(Y[,i],res_alt$U_rrBLUP[,i],use='p'))
      phenix_Y = sapply(1:ncol(Y),function(i) cor(Y[,i],res_alt$Y_phenix[,i],use='p'))
      phenix_U = sapply(1:ncol(Y),function(i) cor(Y[,i],res_alt$U_phenix[,i],use='p'))
      means_U = sapply(1:ncol(Y),function(i) cor(Y[,i],res_mean$U[,i],use='p'))
      results = data.frame(trait=trait.,i=i,SiteYear = colnames(Y),BSFG_Eta,BSFG_U,rrBLUP,phenix_Y,phenix_U,means_U)
    })
  }
}
results$traitID = trait_names[results$trait]

write.csv(results,file = 'Results/collected_results.csv',row.names=F)


library(ggplot2)
library(mapdata)
library(usmap)
library(cowplot)
library(ggthemes)
library(dplyr)
library(tidyr)
source('my_pylo2map.R')
library(geosphere)
library(foreach)
library(data.table)


results = read.csv('Results/collected_results.csv')
traits = unique(results$trait)[c(2,1,3,4)]

summ = aggregate(cbind(BSFG_Eta,BSFG_U,rrBLUP,phenix_Y,phenix_U,means_U)~trait,results,FUN=mean)
summ


results_subset = results[,c('trait','i','SiteYear','BSFG_Eta','phenix_Y','means_U','rrBLUP')]
colnames(results_subset)[4:7] = c('MegaLMM','phenix','GBLUP(env BLUPs)','GBLUP(univariate)')

library(emmeans)
summ_results = foreach(trait. = traits,.combine = 'rbind') %do% {
  results_trait = subset(results_subset,trait == trait.)
  # results_trait_means = aggregate(cbind(BSFG_Eta,BSFG_U,rrBLUP,phenix_Y,phenix_U,means_U)~SiteYear+trait,results_trait,FUN=mean)
  results_trait_means = aggregate(cbind(MegaLMM,phenix,`GBLUP(env BLUPs)`,`GBLUP(univariate)`)~SiteYear+trait,results_trait,FUN=mean)
  results_trait_means_tall = pivot_longer(results_trait_means,cols = -c(1:2))
  res = lm(value ~ SiteYear + name,results_trait_means_tall)
  effects = summary(emmeans(res,spec = 'name'))
  data.frame(Trait = trait.,as.data.frame(effects))
}
summ_results$Trait = factor(summ_results$Trait,levels = unique(summ_results$Trait)[c(2,1,3,4)])
position = position_dodge()
summ_results$name = factor(summ_results$name,levels = unique(summ_results$name)[c(2,1,4,3)])
(p_accuracy <- ggplot(summ_results,aes(x=Trait,y=emmean)) +
    geom_errorbar(aes(group = interaction(name,Trait),ymin = lower.CL,ymax = upper.CL),position=position_dodge(),size=.5) +
    geom_bar(aes(group = interaction(name,Trait),fill = name),stat='identity',position = position_dodge()) +
    scale_fill_manual(name = 'Method',values=scales::brewer_pal(palette = 'Set2')(4)[c(3,4,1,2)]) +
    # scale_color_manual(values = scales::brewer_pal(palette = 'Set2')(4)[c(2,1,3,4)],drop=FALSE)
    ylab('Estimated prediction accuracy') +
    theme_cowplot() + background_grid(major = 'xy') + theme(legend.position = 'bottom') + guides(fill = guide_legend(nrow=2))
)
save_plot('Figures/pred_accuracies.pdf',p_accuracy,base_asp = 1,base_height = 5)


#--------------------------------------#
# Setup map
#--------------------------------------#

fields_metadata = fread('Data/fields_locations_allYears.csv',data.table=F)
fields_metadata$SiteYear = gsub('-','.',fields_metadata$SiteYear)
fields_metadata$Year = sapply(fields_metadata$SiteYear,function(x) {a=strsplit(x,'.',fixed=T)[[1]];a[length(a)]})
fields_metadata$location = paste(fields_metadata$Experiment,fields_metadata$City)
fields_metadata$location[fields_metadata$location == 'ARH1 Mariana'] = 'ARH1 Marianna'

site_distance = distm(fields_metadata[,c('long','lat')])
rownames(site_distance) = colnames(site_distance) = fields_metadata$SiteYear


states = map_data('state')

# adjust map positions
# overlapping locations shifted vertically
# overlapping years shifted horizontally
fields_metadata_start = fields_metadata

experiments = aggregate(cbind(lat,long)~location,data = fields_metadata,FUN = mean)
experiments_start = experiments
dist_y = .5
nudge = T
while(nudge) {
  nudge = F
  experiments_dist = as.matrix(dist(experiments[,c('lat','long')]))
  for(i in 1:nrow(experiments)) {
    j = which(experiments_dist[i,] < dist_y)
    if(length(j) > 1) {
      adj = seq(-1,1,length=length(j))*dist_y*(length(j)-1)/2
      adj = adj[order(experiments[j,]$lat)]
      experiments[j,]$lat = experiments[j,]$lat + adj
      nudge = T
      experiments_dist = as.matrix(dist(experiments[,c('lat','long')]))
    }
  }
}

fields_metadata = fields_metadata_start
fields_metadata$lat = experiments$lat[match(fields_metadata$location,experiments$location)]
fields_metadata$long = experiments$long[match(fields_metadata$location,experiments$location)]
dist_x = .5
for(loc in unique(fields_metadata$location)) {
  j = which(fields_metadata$location == loc)
  if(length(j)>1) {
    adj = seq(-1,1,length=length(j))*dist_x*(length(j)-1)/2
    adj = adj[order(fields_metadata[j,]$long)]
    fields_metadata[j,]$long = fields_metadata[j,]$long + adj
  }
}
(ggplot(states,aes(long,lat)) +
    geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-67))
  # +geom_segment(data = cors_trait2,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor),size=1)
  # +geom_point(data = fields_metadata_start,aes(x=long,y=lat),size=1)
  +geom_point(data = fields_metadata,aes(x=long,y=lat,color = Year),size=1)
)
# (ggplot(states,aes(long,lat)) +
#     geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() 
#     # +coord_cartesian(xlim = c(-104,-67),ylim = c(40,90))#
#   + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-67))
#   # +geom_segment(data = cors_trait2,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor),size=1)
#   # +geom_point(data = fields_metadata_start,aes(x=long,y=lat),size=1)
#   +geom_point(data = fields_metadata,aes(x=long,y=lat),size=2)
#   # +geom_tile(data=na.omit(as.data.frame(r[[1]],xy=T)),aes(fill = bio1,x=x,y=y))
# )
# 
# fields_metadata = fread('Data/fields_locations_allYears.csv',data.table=F)
# states = map_data('state')
# (ggplot(states,aes(long,lat)) +
#     geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() 
#   + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-67))
#   +geom_point(data = fields_metadata,aes(x=long,y=lat),size=2)
# )
# 
# 
# library(raster)
# library(sp)
# 
# r <- getData("worldclim",var="bio",res=10)
# r <- r[[c(1,12)]]
# names(r) <- c("Temp","Prec")
# plot(r[[1]])
# spg = states
# coordinates(spg) <- ~ x + y
# # coerce to SpatialPixelsDataFrame
# gridded(spg) <- TRUE
# # coerce to raster
# rasterDF <- raster(spg)
# rasterDF
# polys = tapply(1:nrow(states),states$group,function(x) Polygons(states[x,1:2],'A'))
# SpatialPolygonsDataFrame(polys)
# r2 = mask(r[[1]],SpatialPolygonsDataFrame(Polygon(states[,1:2])))
# 
# gtf = us[match(toupper(unique(states$region)),toupper(us$NAME_1)),]
# r2 = mask(r[[2]],gtf)
# plot(r2)
          
#--------------------------------------#
# Statistical analysis of max_Gcor
#--------------------------------------#


Gcors = lapply(trait_names,function(trait.) readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.)))
names(Gcors) = traits

library(DescTools)
library(numform)
transform = function(x) FisherZ(x)
# transform = function(x) x
scales::trans_new('fisherz','FisherZ','FisherZInv')
results_sum=foreach(trait. = traits,.combine = 'rbind') %do% {
  results$MvLMM_effect = unlist(results$BSFG_Eta - results$rrBLUP)
  results$MvLMM_effect_vst = unlist(transform(results$BSFG_Eta) - transform(results$rrBLUP))
  results_sum = aggregate(cbind(rrBLUP,BSFG_Eta,MvLMM_effect,MvLMM_effect_vst) ~ SiteYear,subset(results,trait == trait.),FUN = mean)
  results_sum$MvLMM_effect_percent = results_sum$MvLMM_effect/results_sum$rrBLUP

  results_sum$Year = sapply(results_sum$SiteYear,function(x) {a=strsplit(x,'.',fixed=T)[[1]];a[length(a)]})
  Gcor = Gcors[[trait.]]
  Gcor = Gcor[results_sum$SiteYear,results_sum$SiteYear]
  diag(Gcor) = NA
  results_sum$ave_Gcor = rowMeans(Gcor^2,na.rm=T)
  # results_sum$max_Gcor = apply(Gcor^2,1,max,na.rm=T)
  results_sum$max_Gcor = apply((Gcor),1,max,na.rm=T)
  # results_sum$Gcor_90q = apply(Gcor^2,1,function(x) quantile(na.omit(x),.9))
  # h2s = readRDS(sprintf('G2F_byTrait_Gcors/h2s_%s_1.rds',trait.))
  # results_sum$h2 = h2s[results_sum$SiteYear]
  data.frame(Trait = trait.,results_sum)
}
results_sum$Trait = factor(results_sum$Trait,levels = traits)
(p_regression = ggplot(results_sum,aes(x=max_Gcor,y=MvLMM_effect_vst)) +
    geom_point() + geom_smooth(span=1.5) + facet_wrap(~Trait) + #coord_equal() +
    # scale_x_continuous(trans = scales::trans_new('fisherz','FisherZ','FisherZInv'),breaks = scales::extended_breaks(n=10),minor_breaks = NULL,labels = function(x) sprintf('%.1f',x)) +
    scale_x_continuous(limits=c(0.05,.95),trans = scales::trans_new('fisherz','FisherZ','FisherZInv'),breaks = scales::extended_breaks(n=10),minor_breaks = NULL,labels = function(x) f_num(x,digits=1)) +
    xlab('maximum genetic correlation') + #with other site:year
    ylab('Benefit of MvLMM'))  # (difference of z-transformed correlations)
# (p_regression = ggplot(results_sum,aes(x=ave_Gcor,y=MvLMM_effect_vst)) + geom_smooth(span=1.5) + geom_point() + facet_wrap(~Trait) + scale_x_continuous(trans=scales::trans_new('FisherZ',FisherZ,FisherZInv)))
save_plot('Figures/Max_gcov_vs_MvLMM_effect.pdf',p_regression,base_asp = 1,base_height = 8)


models = foreach(trait. = traits) %do% {
  m = lm(MvLMM_effect_vst ~ Year + rrBLUP + max_Gcor,subset(results_sum,Trait == trait. & rrBLUP > -10.2))
  m
}
names(models) = traits
sapply(models,function(x) summary(x)$coef['max_Gcor',])

# plot_grid(p_accuracy+theme(legend.position = 'bottom') + guides(fill = guide_legend(nrow=2)),p_regression,labels = c('A','B'))


#--------------------------------------#
# Map of MvLMM_effect
#--------------------------------------#


limits_Gcor = range(sapply(Gcors,function(x)  range(x[upper.tri(x,diag=F)])))
limits_MvLMM_effect = range(results_sum$MvLMM_effect)
limits_MvLMM_effect = c(0,.4)

MvLMM_maps=foreach(trait. = traits) %do% {
  results_sum_trait = subset(results_sum,Trait == trait.)
  fields_metadata$rrBLUP = results_sum_trait$rrBLUP[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MegaLMM = results_sum_trait$BSFG_Eta[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MvLMM_effect = results_sum_trait$MvLMM_effect[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MvLMM_effect = pmin(limits_MvLMM_effect[2],fields_metadata$MvLMM_effect)
  fields_metadata$MvLMM_effect = pmax(limits_MvLMM_effect[1],fields_metadata$MvLMM_effect)

  (p <- ggplot(states,aes(long,lat)) + ggtitle(trait.) + theme(plot.title = element_text(hjust = 0.5,size = 8)) +
          geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-73))
        # +geom_segment(data = cors_trait2,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor),size=1)
        # +geom_point(data = fields_metadata_start,aes(x=long,y=lat),size=1)
        +geom_point(data = subset(fields_metadata,!is.na(MvLMM_effect)),aes(x=long,y=lat,color = MvLMM_effect,size=MvLMM_effect))
        +scale_size_area(max_size=2)
        + scale_color_viridis_c(name = 'MvLMM\neffect',limits = limits_MvLMM_effect) + guides(size = F) #limits = range(results_sum$MvLMM_effect)
  )
  p

  # Gcor = Gcors[[trait.]]
  # clust = hclust(as.dist(1-Gcor))
  # # print(ggdendro::ggdendrogram(clust) + ggtitle(trait.))
  # library(phytools)
  # coords = fields_metadata[match(rownames(Gcor),fields_metadata$SiteYear),c('lat','long')]
  # coords = as.matrix(coords)
  # colors = fields_metadata$MvLMM_effect[match(rownames(Gcor),fields_metadata$SiteYear)]
  # colors = scales::col_bin(palette = viridis::viridis(100),range(colors))(colors)
  # # colors = factor(sapply(rownames(Gcor),function(x) {a=strsplit(x,'.',fixed=T)[[1]];a[length(a)]}))
  # # rownames(coords) = sapply(rownames(Gcor),function(x) {a=strsplit(x,'.',fixed=T)[[1]];paste(a[-length(a)],sep='_')})
  # rownames(coords) = rownames(Gcor)
  # names(colors) = rownames(coords)
  #
  #
  # a = phylo.to.map(as.phylo.hclust(clust),as.matrix(coords),database='state',plot=F)
  # # a$tree$tip.label = sapply(a$tree$tip.label,function(x) {a=strsplit(x,'.',fixed=T)[[1]];paste(a[-length(a)],sep='_')})
  # my_phylo2map(a,colors = c(colors),split = c(.3,.55),ftype='off',lwd=.8,lty=1,mar=c(0,2,0,0),main = trait.,cex.points=c(.8,1.5)) #,fsize=.4
  #
}


library(DescTools)
library(ggnewscale)
sy_plots=foreach(trait. = traits) %do% {
  print(trait.)
  Gcor = Gcors[[trait.]]
  results_sum_trait = subset(results_sum,Trait == trait.)
  fields_metadata$rrBLUP = results_sum_trait$rrBLUP[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MegaLMM = results_sum_trait$BSFG_Eta[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MvLMM_effect = results_sum_trait$MvLMM_effect[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]

  rank = 1
  foreach(rank = 1) %do% {
    print(rank)
    i = order(-fields_metadata$MvLMM_effect)[rank]
    sy = fields_metadata$SiteYear[i]

    cors = data.frame(SiteYear1 = sy,SiteYear2 = colnames(Gcor),cor = Gcor[sy,])
    cors$Lat1 = fields_metadata$lat[i]
    cors$Lon1 = fields_metadata$long[i]
    cors$Lat2 = fields_metadata$lat[match(colnames(Gcor),fields_metadata$SiteYear)]
    cors$Lon2 = fields_metadata$long[match(colnames(Gcor),fields_metadata$SiteYear)]
    # cors$cor = FisherZ(cors$cor)
    # cors$cor[is.infinite(cors$cor)] = NA#1.1*max(cors$cor[!is.infinite(cors$cor)])
    cors = cors[order(cors$cor),]
    # cors$cor = seq(1:nrow(cors))

    (p <- ggplot(states,aes(long,lat)) + #ggtitle(trait.) +
        geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-73))
      +geom_segment(data = cors,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor))
      # + scale_color_viridis_c(limits = range(cors$cor),name = 'Genetic\ncorrelation')
      + scale_color_viridis_c(limits = limits_Gcor,name = 'Genetic\ncorrelation',option = 'magma')
      # + scale_color_viridis_c(limits = cor_range)
      # +geom_point(data = fields_metadata_start,aes(x=long,y=lat),size=1)
      # + geom_point(data = cors,aes(x=Lon2,y=Lat2,color = cor))
      # +new_scale_color()
      # +geom_point(data = subset(fields_metadata,!is.na(MvLMM_effect)),aes(x=long,y=lat,color = Year))
      # +scale_color_continuous(aesthetics = 'fill')
      # + guides(size = F) #limits = range(results_sum$MvLMM_effect)
      + scale_size(range = c(.1,4)/4,guide=F)
      # + theme(legend.position = 'bottom')
    )
    p
  }
}
sy_plots = do.call(c,sy_plots)

nl = theme(legend.position = 'none')
p = plot_grid(plotlist = lapply(c(MvLMM_maps,sy_plots),function(x) x+nl),ncol=4)
p2 = plot_grid(p,plot_grid(get_legend(MvLMM_maps[[1]]),get_legend(sy_plots[[1]]),ncol=1),rel_widths = c(10,1))
save_plot(filename = 'Figures/MvLMM_effects_top_example.pdf',p2,base_width = 8)


p_dist = foreach(trait. = traits) %do% {
  Gcor = Gcors[[trait.]]
  results_sum_trait = subset(results_sum,Trait == trait.)

  Gcor_tall = data.frame(SY1 = rep(rownames(Gcor),each = ncol(Gcor)),SY2 = rep(colnames(Gcor),nrow(Gcor)),Gcor = c(Gcor))
  dist_tall = data.frame(SY1 = rep(rownames(site_distance),each = ncol(site_distance)),SY2 = rep(colnames(site_distance),nrow(site_distance)),site_distance = c(site_distance)/1e6)
  results = merge(Gcor_tall,dist_tall)
  results = subset(results,SY1 != SY2)
  results$MvLMM_effect = results_sum_trait$MvLMM_effect[match(results$SY1,results_sum_trait$SiteYear)]
  results$MvLMM_effect = pmin(limits_MvLMM_effect[2],results$MvLMM_effect)
  results$MvLMM_effect = pmax(limits_MvLMM_effect[1],results$MvLMM_effect)
  results_sum_trait = results_sum_trait[order(results_sum_trait$MvLMM_effect),]
  results$SY1 = factor(results$SY1,levels = results_sum_trait$SiteYear)

  ggplot(results,aes(x=site_distance,y=Gcor)) + #ylim(limits_Gcor) +
    geom_smooth(aes(group = SY1,color = MvLMM_effect),se=F,span=1,size=.3) +
    scale_color_viridis_c(limits = limits_MvLMM_effect) + xlab('') + ylab('') + theme_cowplot()
}

effect_row = plot_grid(plotlist = lapply(MvLMM_maps,function(x) x+nl),nrow = 1,align = 'v',axis = 'r') + theme(plot.margin = margin(c(0,0,0,20)))
cors_row = plot_grid(plotlist = lapply(sy_plots,function(x) x+nl),nrow = 1) + theme(plot.margin = margin(c(0,0,0,20)))
distances = plot_grid(plotlist = lapply(p_dist,function(x) x+theme(aspect.ratio=1)+nl),nrow = 1)
distances2 = ggdraw(distances) + draw_label('Genetic\ncorrelation',angle = 90,x=-.01,y=.6,vjust = 1.5,size = 10) + draw_label('Distance (1000 Km)',y=0,vjust=-1.5,size = 10) #expression(rho[g])
p1 = plot_grid(effect_row,distances2,cors_row,nrow=3,labels = c('A','B','C'))
lp = theme(legend.justification = 'center',legend.position = 'right')
legends = plot_grid(NULL,get_legend(MvLMM_maps[[1]]+lp),NULL,get_legend(sy_plots[[1]]+lp),NULL,ncol=1,align = 'v',rel_heights = c(1,1,1,1,1))
p2 = plot_grid(p1,legends,nrow=1,rel_widths = c(.9,.1))

save_plot(filename = 'Figures/MvLMM_effects_3_parts_V2.pdf',p2,base_width = 6)


p_dist_mean = foreach(trait. = traits) %do% {
  Gcor = Gcors[[trait.]]
  results_sum_trait = subset(results_sum,Trait == trait.)

  Gcor_tall = data.frame(SY1 = rep(rownames(Gcor),each = ncol(Gcor)),SY2 = rep(colnames(Gcor),nrow(Gcor)),Gcor = c(Gcor))
  dist_tall = data.frame(SY1 = rep(rownames(site_distance),each = ncol(site_distance)),SY2 = rep(colnames(site_distance),nrow(site_distance)),site_distance = c(site_distance)/1e6)
  results = merge(Gcor_tall,dist_tall)
  results = subset(results,SY1 != SY2)
  results$MvLMM_effect = results_sum_trait$MvLMM_effect[match(results$SY1,results_sum_trait$SiteYear)]
  results$MvLMM_effect = pmin(limits_MvLMM_effect[2],results$MvLMM_effect)
  results$MvLMM_effect = pmax(limits_MvLMM_effect[1],results$MvLMM_effect)
  results_sum_trait = results_sum_trait[order(results_sum_trait$MvLMM_effect),]
  results$SY1 = factor(results$SY1,levels = results_sum_trait$SiteYear)

  ggplot(results,aes(x=site_distance,y=Gcor)) + #ylim(limits_Gcor) +
    geom_smooth(se=T,span=1) +
    scale_color_viridis_c(limits = limits_MvLMM_effect) + xlab('') + ylab('') + theme_cowplot()
}
(p3 <- plot_grid(plotlist = lapply(p_dist_mean,function(x) x+theme(aspect.ratio=1)+nl),nrow = 1))
distances = plot_grid(plotlist = lapply(p_dist_mean,function(x) x+theme(aspect.ratio=1)+nl),nrow = 1)
distances2 = ggdraw(distances) + draw_label('Genetic\ncorrelation',angle = 90,x=-.01,y=.6,vjust = 1.5,size = 10) + draw_label('Distance (1000 Km)',y=0,vjust=-1.5,size = 10) #expression(rho[g])
p1 = plot_grid(effect_row,distances2,cors_row,nrow=3,labels = c('A','B','C'))
legends = plot_grid(NULL,get_legend(MvLMM_maps[[1]]+lp),NULL,get_legend(sy_plots[[1]]+lp),NULL,ncol=1,align = 'v',rel_heights = c(1,1,1,1,1))
(p3 <- plot_grid(p1,legends,nrow=1,rel_widths = c(.9,.1)))

save_plot(filename = 'Figures/MvLMM_effects_3_parts_V3.pdf',p3,base_width = 6)


# p_regressions = foreach(trait. = traits) %do% {
#   ggplot(subset(results_sum,Trait == trait.),aes(x=max_Gcor,y=MvLMM_effect_vst)) + ggtitle(trait.) + theme(plot.margin = margin(5,0,0,-20),plot.title = element_text(hjust = 0.5,size = 8)) +
#     # theme()
#     theme(aspect.ratio=1) +
#     geom_point(size=.6) + geom_smooth(span=1.5) + #coord_equal(ylim = range(results_sum$MvLMM_effect_vst)) +
#     scale_x_continuous(limits=c(0.05,.95),trans = scales::trans_new('fisherz','FisherZ','FisherZInv'),breaks = scales::extended_breaks(n=5),minor_breaks = NULL,labels = function(x) f_num(x,digits=1)) +
#     ylim(range(results_sum$MvLMM_effect_vst)) + xlab('') + ylab('')
#     #xlab('maximum genetic correlation') + #with other site:year
#     #ylab('Benefit of MvLMM')
# }
# p_regressions[[1]]
# regressions_row = plot_grid(plotlist = p_regressions,nrow = 1) + theme(plot.margin = margin(c(0,0,10,20))) #+ geom_text(x=0,label = 'adsf',y=10) #+ coord_cartesian(clip = 'off')
# regressions_row2 = ggdraw(regressions_row) + draw_label('Benefit of MvLMM',angle = 90,x=0,vjust = 1.5,size = 8) + draw_label('Maximum genetic correlation',y=0,vjust=-1.5,size = 10)
# effect_row = plot_grid(plotlist = lapply(MvLMM_maps,function(x) x+nl),nrow = 1,align = 'v',axis = 'r') + theme(plot.margin = margin(c(0,0,0,20)))
# cors_row = plot_grid(plotlist = lapply(sy_plots,function(x) x+nl),nrow = 1) + theme(plot.margin = margin(c(0,0,0,20)))
# lp = theme(legend.justification = 'center',legend.position = 'right')
# p3 = plot_grid(regressions_row2,effect_row,cors_row,nrow = 3,rel_heights = c(1.5,1,1),greedy = T,labels = c('A','B','C'))
# lb = get_legend(p_regressions[[1]])
# p4 = plot_grid(p3,plot_grid(lb,get_legend(MvLMM_maps[[1]]+lp),get_legend(sy_plots[[1]]+lp),nrow=3,align = 'v'),ncol=2,rel_widths = c(9,1),greedy = F)
# # p3 = plot_grid(effect_row,cors_row,nrow = 2)
# save_plot(filename = 'Figures/MvLMM_effects_3_parts.pdf',p4,base_width = 5)




op=par(mfrow=c(2,2))
for(trait. in traits[c(1,4,5,8)]) {
  results_trait = subset(results,trait == trait.)
  # plot(results_trait$BSFG_Eta,results_trait$BSFG_U,ylim = c(-.1,1),xlim=c(-.1,1),main = trait.);abline(0,1)
  plot(results_trait$BSFG_Eta,results_trait$means_U,ylim = c(-.1,1),xlim=c(-.1,1));abline(0,1)
  plot(results_trait$BSFG_Eta,results_trait$rrBLUP,ylim = c(-.1,1),xlim=c(-.1,1));abline(0,1)
  plot(results_trait$BSFG_Eta,results_trait$phenix_U,ylim = c(-.1,1),xlim=c(-.1,1));abline(0,1)
  plot(results_trait$BSFG_Eta,results_trait$phenix_Y,ylim = c(-.1,1),xlim=c(-.1,1));abline(0,1)
}
par(op)




aggregate(cbind(BSFG_Eta-rrBLUP,BSFG_U-rrBLUP,phenix_Y-rrBLUP,phenix_U-rrBLUP,means_U-rrBLUP)~trait,results,FUN=mean)
aggregate(cbind(BSFG_Eta-means_U,BSFG_U-means_U,phenix_Y-means_U,phenix_U-means_U,rrBLUP-means_U)~trait,results,FUN=mean)


fields_metadata = fread('fields_locations_allYears.csv',data.table=F)
fields_metadata$SiteYear = gsub('-','.',fields_metadata$SiteYear)
fields_metadata$Year = sapply(fields_metadata$SiteYear,function(x) {a=strsplit(x,'.',fixed=T)[[1]];a[length(a)]})

states = map_data('state')

# adjust map positions
# overlapping locations shifted vertically
# overlapping years shifted horizontally
fields_metadata$location = paste(fields_metadata$Experiment,fields_metadata$City)
fields_metadata$location[fields_metadata$location == 'ARH1 Mariana'] = 'ARH1 Marianna'
fields_metadata_start = fields_metadata

experiments = aggregate(cbind(lat,long)~location,data = fields_metadata,FUN = mean)
experiments_start = experiments
dist_y = .5
nudge = T
while(nudge) {
  nudge = F
  experiments_dist = as.matrix(dist(experiments[,c('lat','long')]))
  for(i in 1:nrow(experiments)) {
    j = which(experiments_dist[i,] < dist_y)
    if(length(j) > 1) {
      adj = seq(-1,1,length=length(j))*dist_y*(length(j)-1)/2
      adj = adj[order(experiments[j,]$lat)]
      experiments[j,]$lat = experiments[j,]$lat + adj
      nudge = T
      experiments_dist = as.matrix(dist(experiments[,c('lat','long')]))
    }
  }
}

# (ggplot(states,aes(long,lat)) +
# geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-67))
# # +geom_segment(data = cors_trait2,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor),size=1)
#   +geom_point(data = experiments_start,aes(x=long,y=lat),size=1)
#   +geom_point(data = experiments,aes(x=long,y=lat),size=1,color='red')
#   )

fields_metadata = fields_metadata_start
fields_metadata$lat = experiments$lat[match(fields_metadata$location,experiments$location)]
fields_metadata$long = experiments$long[match(fields_metadata$location,experiments$location)]
dist_x = .5
for(loc in unique(fields_metadata$location)) {
  j = which(fields_metadata$location == loc)
  if(length(j)>1) {
    adj = seq(-1,1,length=length(j))*dist_x*(length(j)-1)/2
    adj = adj[order(fields_metadata[j,]$long)]
    fields_metadata[j,]$long = fields_metadata[j,]$long + adj
  }
}
(ggplot(states,aes(long,lat)) +
    geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-67))
  # +geom_segment(data = cors_trait2,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor),size=1)
  # +geom_point(data = fields_metadata_start,aes(x=long,y=lat),size=1)
  +geom_point(data = fields_metadata,aes(x=long,y=lat,color = Year),size=1)
)


# trait. = traits[4]
# figs = foreach(trait. = traits[c(1,4,5,8)]) %do% {
#   Gcor = readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.))
#   # print(dim(Gcor))
#   diag(Gcor) = NA
#   Gcor = Gcor^2
#
#
#   fields_metadata$Year = factor(fields_metadata$Year)
#   max_cors = data.frame(SiteYear1 = colnames(Gcor),SiteYear2 = colnames(Gcor)[apply(Gcor,1,which.max)],stringsAsFactors = F)
#   max_cors$lat1 = fields_metadata$lat[match(max_cors$SiteYear1,fields_metadata$SiteYear)]
#   max_cors$long1 = fields_metadata$long[match(max_cors$SiteYear1,fields_metadata$SiteYear)]
#   max_cors$lat2 = fields_metadata$lat[match(max_cors$SiteYear2,fields_metadata$SiteYear)]
#   max_cors$long2 = fields_metadata$long[match(max_cors$SiteYear2,fields_metadata$SiteYear)]
#   max_cors$Year1 = fields_metadata$Year[match(max_cors$SiteYear1,fields_metadata$SiteYear)]
#   max_cors$Year2 = fields_metadata$Year[match(max_cors$SiteYear2,fields_metadata$SiteYear)]
#   # max_cors = data.frame(max_cors,stringsAsFactors = F)
#   all_years <- (ggplot(states,aes(long,lat)) +
#           geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-103,-67))
#         +geom_segment(data = max_cors,aes(x=long1,xend=long2,y=lat1,yend=lat2,color=Year2),size=.5,arrow = arrow(length = unit(0.2,"cm")))
#         +geom_point(data = subset(fields_metadata,SiteYear %in% c(max_cors$SiteYear1,max_cors$SiteYear2)),aes(x=long,y=lat,color = Year),size=1)
#         + scale_color_hue(drop=F,name='Year') + ggtitle(gsub('_',' ',trait.))
#   )# + nl
#   all_years
#
#   # year. = 2015
#   # maps = foreach(year. = 2014:2017) %do% {
#   #   year_cors = subset(max_cors,Year1 == year.)
#   #   year_cors$Y
#   #   p <- (ggplot(states,aes(long,lat)) +
#   #       geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-100,-67))
#   #     +geom_segment(data = year_cors,aes(x=long1,xend=long2,y=lat1,yend=lat2,color=Year2),size=.5,arrow = arrow(length = unit(0.2,"cm")))
#   #     +geom_point(data = subset(fields_metadata,SiteYear %in% c(year_cors$SiteYear1,year_cors$SiteYear2)),aes(x=long,y=lat,color = Year),size=1)
#   #     + scale_color_hue(drop=F) + ggtitle(year.)
#   #   )# + nl
#   #   p
#   # }
#   # # p2 = cowplot::plot_grid(plotlist = maps,nc=2,nr=2)
#   # library(gridExtra)
#   # p2 = grid.arrange(grobs=maps, nrow=2)
#   # save_plot(filename = sprintf('Figures/%s_byYear.pdf',trait.),p2,base_height = 6)
#   all_years
# }
# nl = theme(legend.position = 'none')
# all_traits = plot_grid(
#   plot_grid(
#     plotlist = lapply(figs,function(x) x+nl),ncol=2),
#             get_legend(figs[[1]]+theme(legend.position = 'bottom',
#                                        legend.justification = 'center')),ncol=1,rel_heights = c(.8,.2))
# # all_traits = grid.arrange(grobs=lapply(figs,function(x) x+nl), nrow=2)
# save_plot(filename = sprintf('Figures/allTraits_byYear.pdf',trait.),all_traits,base_height = 8,base_asp = 1)
#
#
# trait. = traits[8]
# maps = foreach(trait. = traits[c(1,4,5,8)]) %do% {
#   Gcor = readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.))
#   diag(Gcor) = NA
#   h2s = readRDS(sprintf('G2F_byTrait_Gcors/h2s_%s_1.rds',trait.))
#   # print(dim(Gcor))
#
#
#   fields_metadata$Year = factor(fields_metadata$Year)
#   fields_metadata$h2 = h2s[fields_metadata$SiteYear]
#   max_cors = data.frame(SiteYear1 = colnames(Gcor),SiteYear2 = colnames(Gcor)[apply(Gcor,1,which.max)],stringsAsFactors = F)
#   max_cors$lat1 = fields_metadata$lat[match(max_cors$SiteYear1,fields_metadata$SiteYear)]
#   max_cors$long1 = fields_metadata$long[match(max_cors$SiteYear1,fields_metadata$SiteYear)]
#   max_cors$lat2 = fields_metadata$lat[match(max_cors$SiteYear2,fields_metadata$SiteYear)]
#   max_cors$long2 = fields_metadata$long[match(max_cors$SiteYear2,fields_metadata$SiteYear)]
#   max_cors$Year1 = fields_metadata$Year[match(max_cors$SiteYear1,fields_metadata$SiteYear)]
#   max_cors$Year2 = fields_metadata$Year[match(max_cors$SiteYear2,fields_metadata$SiteYear)]
#   max_cors$h2s2 = h2s[max_cors$SiteYear2]
#   # max_cors = data.frame(max_cors,stringsAsFactors = F)
#
#   p <- (ggplot(states,aes(long,lat)) +
#           geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-100,-67))
#         +geom_segment(data = max_cors,aes(x=long1,xend=long2,y=lat1,yend=lat2,color=Year2),size=.2,arrow = arrow(length = unit(0.05,"cm")))
#         +geom_point(data = subset(fields_metadata,SiteYear %in% c(max_cors$SiteYear1,max_cors$SiteYear2)),aes(x=long,y=lat,size = h2,color = Year))
#         # +geom_point(data = max_cors,aes(x=long2,y=lat2,color = Year2,size = h2s2))
#         + scale_color_hue(drop=F) + ggtitle(trait.) + scale_radius(limits = c(0,1),range = c(0,1))
#   )# + nl
#   p
# }
# # p2 = cowplot::plot_grid(plotlist = maps,nc=2,nr=2)
# library(gridExtra)
# p2 = grid.arrange(grobs=maps, nrow=2)
# save_plot(filename = ('Figures/byTrait.pdf'),p2,base_height = 6)
#
#
# trait. = traits[c(1,4,5,8)][1]
# Gcor = readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.))
# # print(dim(Gcor))
# diag(Gcor) = NA
#
# SY = unique(fields_metadata$SiteYear)[2]
# pdf(sprintf('Figures/bySite/%s.pdf',trait.))
# foreach(SY = colnames(Gcor) ) %do% {
# cors = data.frame(SiteYear1 = SY,SiteYear2 = colnames(Gcor),cor = Gcor[SY,],stringsAsFactors = F)
# cors$lat1 = fields_metadata$lat[match(cors$SiteYear1,fields_metadata$SiteYear)]
# cors$long1 = fields_metadata$long[match(cors$SiteYear1,fields_metadata$SiteYear)]
# cors$lat2 = fields_metadata$lat[match(cors$SiteYear2,fields_metadata$SiteYear)]
# cors$long2 = fields_metadata$long[match(cors$SiteYear2,fields_metadata$SiteYear)]
# cors$Year1 = fields_metadata$Year[match(cors$SiteYear1,fields_metadata$SiteYear)]
# cors$Year2 = fields_metadata$Year[match(cors$SiteYear2,fields_metadata$SiteYear)]
# (p <- ggplot(states,aes(long,lat)) +
#         geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-100,-67))
#       +geom_segment(data = cors,aes(x=long1,xend=long2,y=lat1,yend=lat2,color=cor),size=.5)
#   # +geom_point(data = subset(fields_metadata,SiteYear %in% colnames(Gcor)),aes(x=long,y=lat),size=1)
#   +geom_point(data = cors,aes(x=long2,y=lat2,fill = cor,color=cor),size=3)
#   + scale_color_gradient2()
#       + ggtitle(SY,trait.)
# )# + nl
# }
# dev.off()


# regions
library(sp)
usmap::.south_region
usmap::.northeast_region
usmap::.midwest_region
usmap::.mountain
usmap::.west_region
fields_coords = SpatialPoints(usmap::usmap_transform(fields_metadata[,c('long','lat')])[,3:4])
fields_metadata$Region = NA
fields_metadata$Region[which(!is.na(over(fields_coords,map2SpatialPolygons(us_map(include = usmap::.south_region),ID=1))))] = 'South'
fields_metadata$Region[which(!is.na(over(fields_coords,map2SpatialPolygons(us_map(include = usmap::.northeast_region),ID=1))))] = 'North East'
fields_metadata$Region[which(!is.na(over(fields_coords,map2SpatialPolygons(us_map(include = usmap::.midwest_region),ID=1))))] = 'Midwest'
fields_metadata$Region[which(!is.na(over(fields_coords,map2SpatialPolygons(us_map(include = usmap::.mountain),ID=1))))] = 'Mountain'
fields_metadata$Region[which(!is.na(over(fields_coords,map2SpatialPolygons(us_map(include = usmap::.west_region),ID=1))))] = 'West'
fields_metadata$Region[which(!is.na(over(fields_coords,map2SpatialPolygons(us_map(include = usmap::.south_region),ID=1))))] = 'Ontario'

fields_metadata$State = substr(fields_metadata$Experiment,1,2)
fields_metadata$State[fields_metadata$State == 'G2'] = 'WI'
fields_metadata$State2 = state.name[match(fields_metadata$State,state.abb)]
fields_metadata$State2[fields_metadata$State == 'ON'] = 'Ontario'
fields_metadata$Region = NA
south = c('AR','GA','SC','TX','NC')
east = c('NY','DE','ON')
midwest = c('OH','MI','IN','IL')
west = c('MO','IA','NE','SD','CO','KS')
north = c('WI','MN')
fields_metadata$Region[fields_metadata$State %in% south] = 'South'
fields_metadata$Region[fields_metadata$State %in% east] = 'North East'
fields_metadata$Region[fields_metadata$State %in% midwest] = 'Midwest'
fields_metadata$Region[fields_metadata$State %in% west] = 'West'
fields_metadata$Region[fields_metadata$State %in% north] = 'North'



# clusters
pdf('Figures/cluster_maps.pdf',width=8,height=5)
for(trait. in traits[c(1,4,5,8)]) {
  Gcor = readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.))
  clust = hclust(as.dist(1-Gcor))
  # print(ggdendro::ggdendrogram(clust) + ggtitle(trait.))
  library(phytools)
  coords = fields_metadata[match(colnames(Gcor),fields_metadata$SiteYear),c('lat','long')]
  coords = as.matrix(coords)
  colors = factor(sapply(rownames(Gcor),function(x) {a=strsplit(x,'.',fixed=T)[[1]];a[length(a)]}))
  # rownames(coords) = sapply(rownames(Gcor),function(x) {a=strsplit(x,'.',fixed=T)[[1]];paste(a[-length(a)],sep='_')})
  rownames(coords) = rownames(Gcor)
  names(colors) = rownames(coords)

  a = phylo.to.map(as.phylo.hclust(clust),as.matrix(coords),database='state',plot=F)
  # a$tree$tip.label = sapply(a$tree$tip.label,function(x) {a=strsplit(x,'.',fixed=T)[[1]];paste(a[-length(a)],sep='_')})
  my_phylo2map(a,colors = c(colors)+1,split = c(.5,.5),ftype='off',lwd=.5,lty=1,mar=c(0,2,0,0),main = trait.,cex.points=c(.8,1)) #,fsize=.4

  colors = viridis(nrow(coords))[rank(coords[,2])]
  names(colors) = rownames(coords)
  data = data.frame(ID = rownames(Gcor),y = coords[,2])
  res = relmatLmer(y~(1|ID),data,relmat = list(ID=Gcor))
  res = as.data.frame(VarCorr(res))$vcov
  my_phylo2map(a,colors = c(colors),split = c(.5,.5),ftype='off',lwd=.5,lty=1,mar=c(0,2,0,0),main = paste(trait.,sprintf('%0.2f',res[1]/sum(res))),cex.points=c(.8,1)) #,fsize=.4

  colors = viridis(nrow(coords))[rank(coords[,1])]
  names(colors) = rownames(coords)
  data = data.frame(ID = rownames(Gcor),y = coords[,1])
  res = relmatLmer(y~(1|ID),data,relmat = list(ID=Gcor))
  res = as.data.frame(VarCorr(res))$vcov
  my_phylo2map(a,colors = c(colors),split = c(.5,.5),ftype='off',lwd=.5,lty=1,mar=c(0,2,0,0),main = paste(trait.,sprintf('%0.2f',res[1]/sum(res))),cex.points=c(.8,1)) #,fsize=.4

  # region
  colors = factor(fields_metadata$Region[match(rownames(Gcor),fields_metadata$SiteYear)])
  names(colors) = rownames(coords)
  my_phylo2map(a,colors = c(colors)+1,split = c(.5,.5),ftype='off',lwd=.5,lty=1,mar=c(0,2,0,0),main = trait.,cex.points=c(.8,1)) #,fsize=.4


}
dev.off()


library(DescTools)
library(psych)
transform = function(x) FisherZ(x)
# transform = function(x) x
results_sum=foreach(trait. = traits[c(1,4,5,8)],.combine = 'rbind') %do% {
  results$MvLMM_effect = unlist(results$BSFG_Eta - results$rrBLUP)
  results$MvLMM_effect_vst = unlist(transform(results$BSFG_Eta) - transform(results$rrBLUP))
  results_sum = aggregate(cbind(rrBLUP,BSFG_Eta,MvLMM_effect,MvLMM_effect_vst) ~ SiteYear,subset(results,trait == trait.),FUN = mean)

  results_sum$Year = sapply(results_sum$SiteYear,function(x) {a=strsplit(x,'.',fixed=T)[[1]];a[length(a)]})
  Gcor = readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.))
  Gcor = Gcor[results_sum$SiteYear,results_sum$SiteYear]
  diag(Gcor) = NA
  results_sum$ave_Gcor = rowMeans(Gcor^2,na.rm=T)
  results_sum$max_Gcor = apply(Gcor^2,1,max,na.rm=T)
  results_sum$Gcor_90q = apply(Gcor^2,1,function(x) quantile(na.omit(x),.9))
  h2s = readRDS(sprintf('G2F_byTrait_Gcors/h2s_%s_1.rds',trait.))
  results_sum$h2 = h2s[results_sum$SiteYear]
  data.frame(Trait = trait.,results_sum)
}
(p = ggplot(results_sum,aes(x=max_Gcor,y=MvLMM_effect_vst)) + geom_point() + geom_smooth(span=1.5) + facet_wrap(~Trait) + coord_equal())
# (p = ggplot(results_sum,aes(x=ave_Gcor,y=MvLMM_effect_vst)) + geom_smooth(span=1.5) + geom_point() + facet_wrap(~Trait) + scale_x_continuous(trans=scales::trans_new('FisherZ',FisherZ,FisherZInv)))
  save_plot('Figures/Max_gcov_vs_MvLMM_effect.pdf',p,base_asp = 1,base_height = 8)


models = foreach(trait. = traits[c(1,4,5,8)]) %do% {
  m = lm(MvLMM_effect_vst ~ Year + rrBLUP + max_Gcor,subset(results_sum,Trait == trait. & rrBLUP > -10.2))
  m
}
names(models) = traits[c(1,4,5,8)]
# sapply(models,function(x) summary(x)$coef['ave_Gcor',])
sapply(models,function(x) summary(x)$coef['max_Gcor',])
# sapply(models,function(x) summary(x)$coef['Gcor_90q',])


  # print(summary(m))
  # print(anova(m))

  # m2 = lm(MvLMM_effect ~ Year,results_sum)
  # results_sum$MvLMM_effect = resid(m2) + coef(m2)[1]
  # op = par()
  # par(mfrow=c(1,2))
  # boxplot(results_sum$MvLMM_effect~results_sum$Year,main = trait.)
  # plot(results_sum$Gcor,results_sum$MvLMM_effect);abline(0,1)
  # # print(cor(results_sum$rrBLUP,results_sum$MvLMM_effect))
  # par(op)

pdf('Figures/MvLMM_effects.pdf')
r=foreach(trait. = traits[c(1,4,5,8)]) %do% {
  results_sum_trait = subset(results_sum,Trait == trait.)
  fields_metadata$rrBLUP = results_sum_trait$rrBLUP[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MegaLMM = results_sum_trait$BSFG_Eta[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MvLMM_effect = results_sum_trait$MvLMM_effect[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]

  print(ggplot(states,aes(long,lat)) + ggtitle(trait.) +
      geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-67))
    # +geom_segment(data = cors_trait2,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor),size=1)
    # +geom_point(data = fields_metadata_start,aes(x=long,y=lat),size=1)
    +geom_point(data = subset(fields_metadata,!is.na(MvLMM_effect)),aes(x=long,y=lat,color = MvLMM_effect,size=MvLMM_effect))
    + scale_color_viridis_c() + guides(size = F) #limits = range(results_sum$MvLMM_effect)
  )

  Gcor = readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.))
  clust = hclust(as.dist(1-Gcor))
  # print(ggdendro::ggdendrogram(clust) + ggtitle(trait.))
  library(phytools)
  coords = fields_metadata[match(rownames(Gcor),fields_metadata$SiteYear),c('lat','long')]
  coords = as.matrix(coords)
  colors = fields_metadata$MvLMM_effect[match(rownames(Gcor),fields_metadata$SiteYear)]
  colors = scales::col_bin(palette = viridis::viridis(100),range(colors))(colors)
  # colors = factor(sapply(rownames(Gcor),function(x) {a=strsplit(x,'.',fixed=T)[[1]];a[length(a)]}))
  # rownames(coords) = sapply(rownames(Gcor),function(x) {a=strsplit(x,'.',fixed=T)[[1]];paste(a[-length(a)],sep='_')})
  rownames(coords) = rownames(Gcor)
  names(colors) = rownames(coords)


  a = phylo.to.map(as.phylo.hclust(clust),as.matrix(coords),database='state',plot=F)
  # a$tree$tip.label = sapply(a$tree$tip.label,function(x) {a=strsplit(x,'.',fixed=T)[[1]];paste(a[-length(a)],sep='_')})
  my_phylo2map(a,colors = c(colors),split = c(.3,.55),ftype='off',lwd=.8,lty=1,mar=c(0,2,0,0),main = trait.,cex.points=c(.8,1.5)) #,fsize=.4

}
dev.off()

library(DescTools)
library(ggnewscale)
traits = traits[c(1,4,5,8)]
Gcors = lapply(traits,function(trait.) readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.)))
names(Gcors) = traits
# pdf('Figures/MvLMM_effects_examples.pdf')
cor_range = range(sapply(Gcors,function(x) x[upper.tri(x)]))
cor_range = c(-1,1)
sy_plots=foreach(trait. = traits) %do% {
  print(trait.)
  Gcor = Gcors[[trait.]]
  results_sum_trait = subset(results_sum,Trait == trait.)
  fields_metadata$rrBLUP = results_sum_trait$rrBLUP[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MegaLMM = results_sum_trait$BSFG_Eta[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MvLMM_effect = results_sum_trait$MvLMM_effect[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]

  rank = 1
  foreach(rank = 1:3) %do% {
    print(rank)
    i = order(-fields_metadata$MvLMM_effect)[rank]
    sy = fields_metadata$SiteYear[i]

    cors = data.frame(SiteYear1 = sy,SiteYear2 = colnames(Gcor),cor = Gcor[sy,])
    cors$Lat1 = fields_metadata$lat[i]
    cors$Lon1 = fields_metadata$long[i]
    cors$Lat2 = fields_metadata$lat[match(colnames(Gcor),fields_metadata$SiteYear)]
    cors$Lon2 = fields_metadata$long[match(colnames(Gcor),fields_metadata$SiteYear)]
    # cors$cor = FisherZ(cors$cor)
    # cors$cor[is.infinite(cors$cor)] = NA#1.1*max(cors$cor[!is.infinite(cors$cor)])
    cors = cors[order(cors$cor),]
    # cors$cor = seq(1:nrow(cors))

    (p <- ggplot(states,aes(long,lat)) + ggtitle(trait.,sy) +
            geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-67))
          +geom_segment(data = cors,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor,size = cor))
          # + scale_color_viridis_c(limits = range(cors$cor))
          + scale_color_viridis_c(limits = cor_range)
          # +geom_point(data = fields_metadata_start,aes(x=long,y=lat),size=1)
          # + geom_point(data = cors,aes(x=Lon2,y=Lat2,color = cor))
          +new_scale_color()
          # +geom_point(data = subset(fields_metadata,!is.na(MvLMM_effect)),aes(x=long,y=lat,color = Year))
          # +scale_color_continuous(aesthetics = 'fill')
          # + guides(size = F) #limits = range(results_sum$MvLMM_effect)
          + scale_size(range = c(.1,4)/4,trans = 'exp')
    )
    p
  }
}
sy_plots = do.call(c,sy_plots)

pdf('Figures/MvLMM_effects_examples.pdf')
nl = theme(legend.position = 'none')
plot_grid(plotlist = lapply(sy_plots[t(matrix(1:length(sy_plots),ncol=4))],function(x) x+nl),ncol=4)
dev.off()


sy_plots=foreach(trait. = traits) %do% {
  print(trait.)
  Gcor = Gcors[[trait.]]
  results_sum_trait = subset(results_sum,Trait == trait.)
  fields_metadata$rrBLUP = results_sum_trait$rrBLUP[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MegaLMM = results_sum_trait$BSFG_Eta[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MvLMM_effect = results_sum_trait$MvLMM_effect[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]

  rank = 1
  foreach(rank = 1) %do% {
    print(rank)
    i = order(-fields_metadata$MvLMM_effect)[rank]
    sy = fields_metadata$SiteYear[i]

    cors = data.frame(SiteYear1 = sy,SiteYear2 = colnames(Gcor),cor = Gcor[sy,])
    cors$Lat1 = fields_metadata$lat[i]
    cors$Lon1 = fields_metadata$long[i]
    cors$Lat2 = fields_metadata$lat[match(colnames(Gcor),fields_metadata$SiteYear)]
    cors$Lon2 = fields_metadata$long[match(colnames(Gcor),fields_metadata$SiteYear)]
    # cors$cor = FisherZ(cors$cor)
    # cors$cor[is.infinite(cors$cor)] = NA#1.1*max(cors$cor[!is.infinite(cors$cor)])
    cors = cors[order(cors$cor),]
    # cors$cor = seq(1:nrow(cors))

    (p <- ggplot(states,aes(long,lat)) + ggtitle(trait.) +
        geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-104,-67))
      +geom_segment(data = cors,aes(x=Lon1,xend=Lon2,y=Lat1,yend=Lat2,color=cor,size = cor))
      # + scale_color_viridis_c(limits = range(cors$cor))
      + scale_color_viridis_c(limits = cor_range)
      # +geom_point(data = fields_metadata_start,aes(x=long,y=lat),size=1)
      # + geom_point(data = cors,aes(x=Lon2,y=Lat2,color = cor))
      +new_scale_color()
      # +geom_point(data = subset(fields_metadata,!is.na(MvLMM_effect)),aes(x=long,y=lat,color = Year))
      # +scale_color_continuous(aesthetics = 'fill')
      # + guides(size = F) #limits = range(results_sum$MvLMM_effect)
      + scale_size(range = c(.1,4)/4,trans = 'log',guide=F)
      + theme(legend.position = 'bottom')
    )
    p
  }
}
sy_plots = do.call(c,sy_plots)
pdf('Figures/MvLMM_effects_top_example.pdf')
nl = theme(legend.position = 'none')
plot_grid(plotlist = lapply(sy_plots,function(x) x+nl),ncol=2)
dev.off()


names(res) = traits[c(1,4,5,8)]
sapply(res,function(x) summary(x)$coef['Gcor',])
sapply(res,function(x) summary(x)$coef['h2',])

results_sum$Year = sapply(results_sum$SiteYear,function(x) {a=strsplit(x,'.',fixed=T)[[1]];a[length(a)]})
results_sum$State = field
results_sum = separate(results_sum,'SiteYear',c('Experiment','Year'))

mymap = readOGR(dsn = 'shape_files/states_21basic/')
pdf('Figures/project_Gcor_taylor.pdf',width = 5)
r=foreach(trait. = traits[c(1,4,5,8)]) %do% {
  results_sum_trait = subset(results_sum,Trait == trait.)
  fields_metadata$rrBLUP = results_sum_trait$rrBLUP[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MegaLMM = results_sum_trait$BSFG_Eta[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]
  fields_metadata$MvLMM_effect = results_sum_trait$MvLMM_effect[match(fields_metadata$SiteYear,results_sum_trait$SiteYear)]

  Gcor = readRDS(sprintf('G2F_byTrait_Gcors/Gcor_%s_1.rds',trait.))
  rownames(coords) == rownames(Gcor)
  input.mat = Gcor
  pts = SpatialPointsDataFrame(coords[,c('long','lat')],data = data.frame(coords))
  library(recluster)


  #DISTANCE MATRIX
  dis.mat <- as.dist(input.mat)

  #RUN NMDS
  nmds <- metaMDS(dis.mat, distance = "euclidean", trymax=10)

  # Assign a colour to each site according to its position in the ordination space.
  ordi.col<-recluster.col(nmds$points)
  plot.col <- rgb(ordi.col[,c(3,4,5)],maxColorValue=255)

  # par(mfrow=c(2,1), mai=c(.1,.1,.5,.1))
  #PLOT NMDS RESULTS
  plot(nmds$points, col=NULL,bg=plot.col, pch=21, cex=1.5, xlab="NMDS 1", ylab="NMDS 2",main = trait.)
  #THE STRESS VALUE IS HOW WELL THE NMDS RECAPITULATES THE ORIGINAL DISTANCE MATRIX
  #legend("bottomright", legend=paste("stress:", round(nmds_edat$stress, 3), sep=" "))

  #PLOT MAP WITH MATCHING COLOR
  plot(mymap, col=grey(.9,.3),xlim = c(-105,-60),ylim = c(10,50))
  # plot(pts, col=NULL,bg=plot.col, pch=22,cex=2, add=T)
  plot(pts, col=NULL,bg=plot.col, pch=21,cex=1, add=T)
  # plot(pts, col=plot.col, pch=15,cex=2, add=T)


}
dev.off()



# Lambdas
library(BSFG)
trait. = traits[4]
pdf('Figures/Factor_loadings.pdf')
foreach(trait. = traits[c(1,4,5,8)]) %do% {
  Lambda = readRDS(sprintf('G2F_byTrait_Gcors/Lambda_%s_1.rds',trait.))
  fields_metadata_trait = fields_metadata[match(rownames(Lambda),fields_metadata$SiteYear),]
  Lambda = sweep(Lambda,2,sign(colMeans(Lambda)),'*')
  for(i in 1:6) {
    # Image(Lambda[,1:10])
    fields_metadata_trait$lambda = Lambda[,i]
    all_years <- (ggplot(states,aes(long,lat)) +
                    geom_polygon(aes(group=group),fill=NA,col=1,size=.2) + theme_map() + coord_map("albers",  lat0 = 90, lat1 = 40,xlim = c(-103,-67))
                  # +geom_segment(data = fields_metadata_trait,aes(x=long1,xend=long2,y=lat1,yend=lat2,color=Year2),size=.5,arrow = arrow(length = unit(0.2,"cm")))
                  +geom_point(data = fields_metadata_trait,aes(x=long,y=lat,color = lambda),size=3)
                  +scale_color_gradient2(mid='grey70')
                  # + scale_color_hue(drop=F,name='Year')
                  + ggtitle(gsub('_',' ',trait.),i)
    )# + nl
    print(all_years)
  }
}
dev.off()
