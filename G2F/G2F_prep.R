library(data.table)
library(lme4)
library(tidyr)

# Select one GBS sample per inbred line
geno_info = fread('Data/g2f_2014_gbs_data.csv',data.table=F)
taxa_list = fread('Data/g2f_2014_zeaGBSv27_TaxaList.txt',data.table=F)
taxa_summary = fread('Data/g2f_2014_zeaGBSv27_TaxaSummary.txt',data.table=F)
taxa_list$PropMissing = taxa_summary$`Proportion Missing`[match(taxa_list$Taxa,taxa_summary$`Taxa Name`)]
taxa_list$DNAName = geno_info$`<DNAName>`[match(taxa_list$Taxa,geno_info$`<GBSSampleName>`)]

taxa_list$DNAName[taxa_list$DNAName=='Tx303'] = 'TX303'
taxa_list$DNASample[taxa_list$DNASample=='B73_PHG39-15'] = 'B73xPHG39-15'
taxa_list$DNASample[taxa_list$DNASample=='PHN11_OH43_0036'] = 'PHN11_Oh43_0036'
taxa_list$DNASample[taxa_list$DNASample=='PH207_PHG47-5'] = 'PH207xPHG47-5'
taxa_list$DNASample[taxa_list$DNASample=='NYH-258'] = 'NyH-258'
taxa_list$DNASample[taxa_list$DNASample=='MO44_PHW65_0035'] = 'Mo44_PHW65_0035'
taxa_list$DNASample[taxa_list$DNASample=='MO17'] = 'Mo17'
taxa_list$DNASample[taxa_list$DNASample=='554353-1-1-B'] = 'BSSSC0_044'

taxa_list_selected = c()
DNASamples = unique(taxa_list$DNASample)
for(dna in DNASamples) {
  i = taxa_list$DNASample == dna
  taxa_list_i = subset(taxa_list,i)
  taxa_list_i = taxa_list_i[order(taxa_list_i$PropMissing)[1],]
  taxa_list_selected = rbind(taxa_list_selected,tibble(DNASample = dna,DNAName = taxa_list_i$DNAName[1],Taxa = taxa_list_i$Taxa[1]))
}

write.csv(taxa_list_selected,file = 'Data/g2f_2014_zeaGBSv27_DNAtoTaxa.txt',row.names=F)


# load phenotype data
pheno_data_2014 = fread('Data/g2f_2014_hybrid_no_outliers.csv',data.table=F)
pheno_data_2015 = fread('Data/g2f_2015_hybrid_data_clean.csv',data.table=F)
pheno_data_2016 = fread('Data/g2f_2016_hybrid_data_clean.csv',data.table=F)
pheno_data_2017 = fread('Data/g2f_2017_hybrid_data_clean.csv',data.table=F)


# adjust column names of 2014 data
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Rep'] ="Replicate"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Plot area'] ="Plot Area"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Plant height [cm]'] ="Plant Height [cm]"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Ear height [cm]'] = "Ear Height [cm]"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Root lodging [plants]'] = "Root Lodging [plants]"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Stalk lodging [plants]'] = "Stalk Lodging [plants]"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Grain Moisture [percent]'] = "Grain Moisture [%]"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Grain yield [bu/A]'] = "Grain Yield [bu/A]"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Plot Discarded'] = "Plot Discarded [enter \"\"yes\"\" or \"\"blank\"\"]"
colnames(pheno_data_2014)[colnames(pheno_data_2014) == 'Filler'] = "Filler [enter \"\"filler\"\" or \"\"blank\"\"]"
pheno_data_2014$`LOCAL_CHECK (Yes, No[Blank])` = ifelse(pheno_data_2014$`Heterotic Pool` == 0,'Yes','')
pheno_data_2014$Year = 2014


for(col in colnames(pheno_data_2015)[colnames(pheno_data_2015) %in% colnames(pheno_data_2014) == F]) pheno_data_2014[[col]] = NA
pheno_data = rbind(pheno_data_2014[,colnames(pheno_data_2015)],
                   pheno_data_2015,pheno_data_2016,pheno_data_2017)



# find DNAName
pheno_data = separate(pheno_data,'Pedigree',into = c('DNAName','Tester'),sep='/',remove=FALSE)
pheno_data$DNAName_mod = pheno_data$DNAName
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX740'] = "TX740-B5"  # not 100% sure
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX772W'] = "TX772W-B4-B8"  # not 100% sure
pheno_data$DNAName_mod[pheno_data$DNAName_mod == '(TX739);LAMA2002-10-1-B-B-B-B3-B7ORANGE-B6'] = "(TX739)_LAMA2002-10-1-B-B-B-B3-B7_ORANGE-B" # (TX739)LAMA2002-10-1-B-B-B-B3-B7ORANGE-B7-B11 is other option, and very closely related
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX736((TX772XT246)XTX772)-1-5....'] = "(TX736)_((TX772_X_T246)_X_TX772)-1-5-B-B-B-B-B-B6-B6-B2-B13"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == '???TX205'] = "TX205"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == '(LAMA2002-35-2-B-B-B-B_CG44)-1-3-B-1-1-B24-B5-B16'] = "(LAMA2002-35-2-B-B-B-B/CG44)-1-3-B-1-1-B24-B5-B16-B18-B23"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'LAMA2002-12-1-B-B-B-B-B-B-1-B28-B14'] = "LAMA2002-12-1-B-B-B-B-B-B-1-B28-B14-B27-B"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX130 RP'] = "TX130 RP" # unclear if TX130RP-B27-B18 or TX130RP-B15-B5
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX745 RP'] = "TX745RP-B3-B11"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX803'] = "TX745RP-B3-B11"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX768'] = "TX768-B-B18-B18"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX770'] = "TX770RP-B25"  # not 100% sure also TX770RP-B25-B25, but very similar
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX601'] = "TX601-B-B20-B13"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == '(CML442-B*CML343-B-B-B-B-B-B)-B-B-1-1-B-B-B-1-B12-1-B19'] = "(CML442-B/CML343-B-B-B-B-B-B)-B-B-1-1-B-B-B-1-B12-1-B21" # very similar
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'CML450*TX110'] = "CML450-B/TX110)-B-3-B-1-B-B-1-1-B18-B21-B8-1-B25" # very similar
pheno_data$DNAName_mod[grep('MBNIL B',pheno_data$DNAName_mod)] = sub('MBNIL B','MBNILB',pheno_data$DNAName_mod[grep('MBNIL B',pheno_data$DNAName_mod)])
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TR 9-1-1-6'] = "TR9-1-1-6"
pheno_data$DNAName_mod[pheno_data$DNAName_mod == 'TX6252'] = "TX6252-B23"  # not 100% sure
pheno_data$DNAName_mod[pheno_data$DNAName_mod == '78010'] = "DK78010"

pheno_data$DNAName = pheno_data$DNAName_mod

pheno_data$Heterotic_Pool = pheno_data_2014$`Heterotic Pool`[match(pheno_data$Pedigree,pheno_data_2014$Pedigree)]

# remove bad data
pheno_data = subset(pheno_data, `Plot Discarded [enter ""yes"" or ""blank""]` != 'Yes')
pheno_data = subset(pheno_data,Comments == '')

# add ASI
pheno_data$`ASI [days]` = pheno_data$`Pollen DAP [days]` - pheno_data$`Silk DAP [days]`


write.csv(pheno_data,file = 'Data/g2f_2014_2017_hybrid_data_clean.csv',row.names=F)



# field metadata
fields_2014 = fread('Data/g2f_2014_field_characteristics.csv',data.table=F)[,c('Experiment','lat','long','Type','City')]
fields_2015 = fread('Data/g2f_2015_field_metadata.csv',data.table=F)[,c('Experiment','corner1 lat','corner1 lon','Type','City')]
fields_2016 = fread('Data/g2f_2016_field_metadata.csv',data.table=F)[,c('Experiment_Code','Latitude_of_Field_Corner_#1 (lower left)','Longitude_of_Field_Corner_#1 (lower left)','Treatment','City')]
fields_2017 = fread('Data/g2f_2017_field_metadata.csv',data.table=F)[,c('Experiment_Code','Latitude_of_Field_Corner_#1 (lower left)','Longitude_of_Field_Corner_#1 (lower left)','Treatment','City')]
colnames(fields_2015) = c('Experiment','lat','long','Type','City')
colnames(fields_2016) = c('Experiment','lat','long','Type','City')
colnames(fields_2017) = c('Experiment','lat','long','Type','City')


subset(fields_2017,Experiment == 'TXH1')
fields_2017[fields_2017$Experiment == 'TXH1',]$lat = 30.54684
fields_2017[fields_2017$Experiment == 'TXH1',]$long = -96.43468

fields_metadata = rbind(data.frame(SiteYear = paste0(fields_2014$Experiment,'-2014'),fields_2014),
                        data.frame(SiteYear = paste0(fields_2015$Experiment,'-2015'),fields_2015),
                        data.frame(SiteYear = paste0(fields_2016$Experiment,'-2016'),fields_2016),
                        data.frame(SiteYear = paste0(fields_2017$Experiment,'-2017'),fields_2017))
fields_metadata = subset(fields_metadata,toupper(Type) != 'INBRED')



subset(fields_metadata,Experiment == 'IAH2')
fields_metadata[fields_metadata$Experiment == 'IAH2' & is.na(fields_metadata$lat),]$lat = 42.06621
fields_metadata[fields_metadata$Experiment == 'IAH2' & is.na(fields_metadata$long),]$long = -94.72761
subset(fields_metadata,Experiment == 'NYH1')
fields_metadata[fields_metadata$Experiment == 'NYH1' & is.na(fields_metadata$lat),]$lat = 42.72877
fields_metadata[fields_metadata$Experiment == 'NYH1' & is.na(fields_metadata$long),]$long = -76.65165
subset(fields_metadata,Experiment == 'TXH2')
fields_metadata[fields_metadata$Experiment == 'TXH2' & is.na(fields_metadata$lat),]$lat = 34.18467
fields_metadata[fields_metadata$Experiment == 'TXH2' & is.na(fields_metadata$long),]$long = -101.9494
subset(fields_metadata,Experiment == 'IAH1')
fields_metadata[fields_metadata$Experiment == 'IAH1' & is.na(fields_metadata$lat),]$lat = 41.19943
fields_metadata[fields_metadata$Experiment == 'IAH1' & is.na(fields_metadata$long),]$long = -91.49516
subset(fields_metadata,Experiment == 'IAH3')
fields_metadata[fields_metadata$Experiment == 'IAH3' & is.na(fields_metadata$lat),]$lat = 41.97575
fields_metadata[fields_metadata$Experiment == 'IAH3' & is.na(fields_metadata$long),]$long = -92.24144
subset(fields_metadata,Experiment == 'IAH4')
fields_metadata[fields_metadata$Experiment == 'IAH4' & is.na(fields_metadata$lat),]$lat = 41.99439
fields_metadata[fields_metadata$Experiment == 'IAH4' & is.na(fields_metadata$long),]$long = -93.68863
subset(fields_metadata,Experiment == 'NEH1')
fields_metadata[fields_metadata$Experiment == 'NEH1' & is.na(fields_metadata$lat),]$lat = 41.16200
fields_metadata[fields_metadata$Experiment == 'NEH1' & is.na(fields_metadata$long),]$long = -96.40900
subset(fields_metadata,Experiment == 'NEH4')
fields_metadata[fields_metadata$Experiment == 'NEH4' & is.na(fields_metadata$lat),]$lat = 41.16200
fields_metadata[fields_metadata$Experiment == 'NEH4' & is.na(fields_metadata$long),]$long = -96.40900
subset(fields_metadata,Experiment == 'SCH1')
fields_metadata[fields_metadata$Experiment == 'SCH1' & is.na(fields_metadata$lat),]$lat = 34.62261
fields_metadata[fields_metadata$Experiment == 'SCH1' & is.na(fields_metadata$long),]$long = -82.73796
subset(fields_metadata,Experiment == 'ILH1')
fields_metadata[fields_metadata$Experiment == 'ILH1' & is.na(fields_metadata$lat),]$lat = 40.06119
fields_metadata[fields_metadata$Experiment == 'ILH1' & is.na(fields_metadata$long),]$long = -88.23327
subset(fields_metadata,Experiment == 'ILH2')
fields_metadata[fields_metadata$Experiment == 'ILH2' & is.na(fields_metadata$lat),]$lat = 40.08489
fields_metadata[fields_metadata$Experiment == 'ILH2' & is.na(fields_metadata$long),]$long = -88.22541
subset(fields_metadata,Experiment == 'INH1')
fields_metadata[fields_metadata$Experiment == 'INH1' & is.na(fields_metadata$lat),]$lat = 40.47835
fields_metadata[fields_metadata$Experiment == 'INH1' & is.na(fields_metadata$long),]$long = -86.99013

subset(fields_metadata,is.na(lat))

fields_metadata = rbind(fields_metadata,data.frame(SiteYear = c('IAH1a-2014','IAH1b-2014','IAH1c-2014'),
                                                   Experiment = c('IAH1a','IAH1b','IAH1c'),
                                                   subset(fields_metadata,SiteYear == 'IAH1-2014')[,-c(1:2)]))
fields_metadata = rbind(fields_metadata,data.frame(SiteYear = c('WIH1.2014'),
                                                   Experiment = c('WIH1'),
                                                   subset(fields_metadata,SiteYear == 'WIH1-2015')[,-c(1:2)]))
fields_metadata = rbind(fields_metadata,data.frame(SiteYear = c('TXH1-Dry-2017','TXH1-Early-2017','TXH1-Late-2017'),
                                                   Experiment = c('TXH1-Dry','TXH1-Early','TXH1-Late'),
                                                   subset(fields_metadata,SiteYear == 'TXH1-2017')[,-c(1:2)]))

write.csv(fields_metadata,file = 'Data/fields_locations_allYears.csv',row.names=F)

