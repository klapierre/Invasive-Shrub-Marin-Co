library(ggplot2)
library(grid)
library(tidyr)
library(dplyr)

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


setwd('C:\\Users\\Kim\\Dropbox\\2015_MarinParks_LaPierre\\marin parks data')

trt <- read.csv('La Pierre_MPN_treatments_2015.csv')
pots <- read.csv('La Pierre_MPN_pots_2015.csv')
nods <- read.csv('La Pierre_MPN_harvest_2015.csv')%>%
  filter(notes!='no growth', notes!='dead')

nods <- nods%>%
  #combine nodule number and treatment data
  left_join(pots)%>%
  left_join(trt)%>%
  mutate(func_nods=as.numeric(func_nods), nonf_nods=as.numeric(nonf_nods))%>%
  #calculate total nodules
  mutate(total_nods=func_nods+nonf_nods)%>%
  #drop irrelevent information for now
  select(pot, func_nods, nonf_nods, total_nods, harvest.date, planting.date, species, dilution, soil_site, soil_spp, soil_trt)%>%
  filter(dilution!='ctl')


###calculating Most Probable Number (MPN)

#table of number of pots nodulated and number of pots possible to be nodulated
nodsBinary <- nods%>%
  select(-func_nods, -nonf_nods)%>%
  #create binary nodulation state variable
  mutate(nodulated=as.numeric(ifelse(total_nods==0, 0, 1)))%>%
  select(-total_nods)%>%
  #calculate number of plants nodulated and total plants possible (N)
  #note, a few missing individuals from GEMO, LUNA, and ACGL at a few sites bring down their reps from 5 to 4.5 or 4.75 across the dilutions, see below. this is a big problem for LUBI, which should maybe be dropped
  group_by(species, soil_site, soil_spp, soil_trt, dilution)%>%
  summarise(nodulated=sum(nodulated), N=n())

#table of number of nodulated pots for each dilution (in columns)
nodsDilution <- nodsBinary%>%
  select(-N)%>%
  #transpose to spread dilutions
  spread(key=dilution, value=nodulated)

#use BAM Appendix 2 MPN calculator (excel workbook)
MPN <- read.csv('La Pierre_MPN_MPN estimates_2015.csv')

###TO DO
# add stats, standardize variables to uninvaded as controls


###figures of invaded vs uninvaded areas

#GEMO vs uninvaded
gemoMPN <- MPN%>%
  filter(soil_spp=='GEMO' | soil_trt=='uninvaded')%>%
  filter(species!='CYSC', species!='SPJU')%>%
  mutate(soil_trt=factor(soil_trt, levels=c('uninvaded', 'untreated', 'herbicided', 'mowed', 'pulled')))

#CYSC vs uninvaded
cyscMPN <- MPN%>%
  filter(soil_spp=='CYSC' | soil_trt=='uninvaded')%>%
  filter(species!='GEMO', species!='SPJU') 

#SPJU vs uninvaded
spjuMPN <- MPN%>%
  filter(soil_spp=='SPJU' | soil_trt=='uninvaded')%>%
  filter(species!='GEMO', species!='CYSC') 

gemoInvPlot <- ggplot(data=barGraphStats(data=subset(gemoMPN, soil_trt=='uninvaded' | soil_trt=='untreated'), variable="MPN", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('Most Probable Number (MPN)') +
  theme(legend.position=c(0.8,0.9)) +
  scale_y_continuous(limits=c(0,2200)) +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_discrete(labels=c('uninvaded', 'invaded'))

cyscInvPlot <- ggplot(data=barGraphStats(data=cyscMPN, variable="MPN", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('') +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,2200)) +
  scale_x_discrete(limits=c('CYSC', 'ACGL', 'LUBI', 'LUNA'))

spjuInvPlot <- ggplot(data=barGraphStats(data=spjuMPN, variable="MPN", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('') +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,2200)) +
  scale_x_discrete(limits=c('SPJU', 'ACGL', 'LUBI', 'LUNA'))

#3 panel figure
pushViewport(viewport(layout=grid.layout(1,3)))
print(gemoInvPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(cyscInvPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(spjuInvPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3)) 



#species responses to GEMO invasions and treatments
ggplot(data=barGraphStats(data=gemoMPN, variable="MPN", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('Most Probable Number (MPN)') +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  theme(legend.position=c(0.9,0.9))










###looking at number of nodules in highest soil dilution

#subset out the 10^-2 dilution
dil0 <- nods%>%
  filter(dilution==0)

###subset by species
#this is not quite right, because it includes sites where each species isn't found
gemoDil0 <- dil0%>%
  filter(soil_spp=='GEMO' | soil_trt=='uninvaded')%>%
  filter(species!='CYSC', species!='SPJU')%>%
  mutate(soil_trt=factor(soil_trt, levels=c('uninvaded', 'untreated', 'herbicided', 'mowed', 'pulled')))

cyscDil0 <- dil0%>%
  filter(soil_spp=='CYSC' | soil_trt=='uninvaded')%>%
  filter(species!='GEMO', species!='SPJU') 

spjuDil0 <- dil0%>%
  filter(soil_spp=='SPJU' | soil_trt=='uninvaded')%>%
  filter(species!='GEMO', species!='CYSC') 

###quick plots

#species responses to GEMO invasions
ggplot(data=barGraphStats(data=subset(dil0, soil_trt=='uninvaded' | soil_trt=='untreated'), variable="total_nods", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('Total Nodules') +
  theme(legend.position=c(0.8,0.93)) +
  scale_x_discrete(limits=c('GEMO', 'CYSC', 'SPJU', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_discrete(labels=c('uninvaded', 'invaded'))  
  
# gemoInvPlot <- ggplot(data=barGraphStats(data=subset(gemoDil0, soil_trt=='uninvaded' | soil_trt=='untreated'), variable="total_nods", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
#   xlab('Plant Species') +
#   ylab('Total Nodules') +
#   theme(legend.position=c(0.8,0.93)) +
#   scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
#   scale_fill_discrete(labels=c('uninvaded', 'invaded'))
# 
# #species responses to CYSC invasions
# cyscInvPlot <- ggplot(data=barGraphStats(data=cyscDil0, variable="total_nods", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
#   xlab('Plant Species') +
#   ylab('Total Nodules') +
#   theme(legend.position='none') +
#   scale_y_continuous(limits=c(0,20)) +
#   scale_x_discrete(limits=c('CYSC', 'ACGL', 'LUBI', 'LUNA'))
# 
# #species responses to SPJU invasions
# spjuInvPlot <- ggplot(data=barGraphStats(data=spjuDil0, variable="total_nods", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
#   xlab('Plant Species') +
#   ylab('Total Nodules') +
#   theme(legend.position='none') +
#   scale_y_continuous(limits=c(0,20)) +
#   scale_x_discrete(limits=c('SPJU', 'ACGL', 'LUBI', 'LUNA'))
# 
# 
# #3 panel figure
# pushViewport(viewport(layout=grid.layout(1,3)))
# print(gemoInvPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(cyscInvPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(spjuInvPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3))


#species responses to GEMO invasions and treatments
ggplot(data=barGraphStats(data=gemoDil0, variable="total_nods", byFactorNames=c("species", "soil_trt")), aes(x=species, y=mean, fill=soil_trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('Total Nodules') +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  theme(legend.position=c(0.9,0.9))

