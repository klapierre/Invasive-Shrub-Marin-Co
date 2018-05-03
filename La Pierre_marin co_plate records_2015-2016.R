library(ggplot2)
library(grid)
library(tidyverse)
library(lme4)


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}


###########################################################################
###########################################################################

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\Marin County Invasion\\marin parks data')

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\Marin County Invasion\\marin parks data')

#REP PCR record
trt <- read.csv('La Pierre_MPN_treatments_2015.csv')
wholesoilPots <- read.csv('La Pierre_whole soil_pots_2016.csv')%>%
  mutate(experiment='whole_soil', dilution='NA')
MPNpots <- read.csv('La Pierre_MPN_pots_2015.csv')%>%
  mutate(experiment='MPN', planting_date='NA')%>%
  select(rack, pot, species, planting_date, experiment, dilution)
repMPN <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\Marin County Invasion\\marin parks data\\DNA sequences\\PCR\\REP_PCR_strains completed.csv')%>%
  filter(experiment=='MPN')%>%
  left_join(MPNpots)
repWholeSoil <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\Marin County Invasion\\marin parks data\\DNA sequences\\PCR\\REP_PCR_strains completed.csv')%>%
  filter(experiment=='whole_soil')%>%
  left_join(wholesoilPots)

rep <- rbind(repMPN, repWholeSoil)%>%
  left_join(trt)
#write.csv(rep, 'La Pierre_2018_REP isolates.csv')






# #if re-sequencing
# ###data
# trt <- read.csv('La Pierre_MPN_treatments_2015.csv')
# wholesoilPots <- read.csv('La Pierre_whole soil_pots_2016.csv')%>%
#   mutate(experiment='whole_soil', dilution='NA')
# wholesoilNodules <- read.csv('La Pierre_2017_MC_whole soil_plate record_for analysis.csv')%>%
#   mutate(pot=plant_number, freezer_stock=freezer.stock)%>%
#   select(-plant_number, -freezer.stock, -initial_notes, -R1_notes, -R2_notes, -liquid_notes)%>%
#   left_join(wholesoilPots)
# MPNpots <- read.csv('La Pierre_MPN_pots_2015.csv')%>%
#   mutate(experiment='MPN', planting_date='planting.date')%>%
#   select(rack, pot, species, planting_date, experiment, dilution)
# MPNnodules <- read.csv('La Pierre_2016_MC_MPN_plate record_for analysis.csv')%>%
#   mutate(pot=plant_number, R1_plate=R1, R2_plate=R2)%>%
#   select(-plant_number, -notes_R1, -notes_R2, -notes_liquid, -notes_initial, -R1, -R2)%>%
#   left_join(MPNpots)
# 
# nodules <- rbind(wholesoilNodules, MPNnodules)%>%
#   left_join(trt)
# 
# toSequence <- nodules%>%
#   #only keep those that are from uninvaded or untreated soils
#   filter(soil_trt %in% c('uninvaded', 'untreated'))%>%
#   #filter out those we don't have a culture in freezer stock
#   filter(freezer_stock!=0)%>%
#   group_by(species, soil_site, soil_spp, soil_trt)%>% 
#   mutate(random_choice=sample(1:length(species), n(), replace=F))%>%
#   ungroup()%>%
#   mutate(isolate_number=paste(pot, nodule_rep, sep=''))%>%
#   select(isolate_number, pot, nodule_rep, experiment, species, soil_site, soil_spp, soil_trt, random_choice)
# #write.csv(toSequence, 'La Pierre_Marin Co_isolates to sequence.csv', row.names=F)  
# 
# summaryParkToSequence <- toSequence%>%
#   group_by(species, soil_site, soil_spp, soil_trt)%>%
#   summarize(count=length(pot))%>%
#   ungroup()
# 
# summarySpeciesToSequence <- toSequence%>%
#   group_by(species, soil_spp, soil_trt)%>%
#   summarize(count=length(pot))%>%
#   ungroup()









