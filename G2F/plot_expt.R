library(ggplot2)

files = list.files(path='BLUP_matrices/',pattern = 'csv')
traits = sub('_BLUPs.csv','',files,fixed=T)
traits = gsub(' ','_',traits)
traits = traits[c(4,1,2,3)]
trait_names = traits
traits = c('DTS','ASI','Grain Yield','Plant Height')
names(trait_names) = traits