library(ggplot2)
library(grid)
library(tidyr)
library(dplyr)
# library(sjPlot) #to plot lmer effects, doesn't work for some reason
# library(sjmisc) #to plot lmer effects, doesn't work for some reason
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
MPN <- read.csv('La Pierre_MPN_MPN estimates_2015.csv') %>%
  mutate(ln_MPN=log(MPN+1))


###MPN difference
MPNinv <- MPN%>%
  mutate(soil_trt_spp=paste(soil_trt, soil_spp, sep='_'))%>%
  select(-MPN, -MPN_lower, -MPN_upper, -soil_spp, -soil_trt)%>%
  #compare by trt
  spread(key=soil_trt_spp, value=ln_MPN)%>%
  #relativize to uninvaded
  mutate(cysc_effect=ifelse(untreated_CYSC=='NA', 'NA', (untreated_CYSC - uninvaded_uninvaded)))%>%
  mutate(spju_effect=ifelse(untreated_SPJU=='NA', 'NA', (untreated_SPJU - uninvaded_uninvaded)))%>%
  mutate(gemo_effect=ifelse(untreated_GEMO=='NA', 'NA', (untreated_GEMO - uninvaded_uninvaded)))%>%
  mutate(herb_effect=ifelse(herbicided_GEMO=='NA', 'NA', (herbicided_GEMO - uninvaded_uninvaded)))%>%
  mutate(mow_effect=ifelse(mowed_GEMO=='NA', 'NA', (mowed_GEMO - uninvaded_uninvaded)))%>%
  mutate(pull_effect=ifelse(pulled_GEMO=='NA', 'NA', (pulled_GEMO - uninvaded_uninvaded)))

#long form
MPNdiff <- MPNinv%>%
  select(-untreated_CYSC, -untreated_SPJU, -untreated_GEMO, -uninvaded_uninvaded, -herbicided_GEMO, -mowed_GEMO, -pulled_GEMO)%>%
  gather(key=soil_trt_spp, value=MPN_diff, cysc_effect:pull_effect, na.rm=T)


###mixed effects models
#gemo invasion
MPNgemo <- MPNinv%>%
  select(species, soil_site, herbicided_GEMO, mowed_GEMO, pulled_GEMO, uninvaded_uninvaded, untreated_GEMO)%>%
  filter(species!='SPJU', species!='CYSC')%>%
  gather(key=soil_trt_spp, value=ln_MPN, herbicided_GEMO:untreated_GEMO, na.rm=T)

gemoModel <- lmer(ln_MPN ~ species*soil_trt_spp + (1 | soil_site), data=subset(MPNgemo, soil_trt_spp=='uninvaded_uninvaded'|soil_trt_spp=='untreated_GEMO'))
summary(gemoModel)
ranef(gemoModel)
fixef(gemoModel)
# sjp.lmer(gemoModel, type='fe') #doesn't work for some reason

summary(gemoModelglm <- glm(ln_MPN ~ species*soil_trt_spp, data=subset(MPNgemo, soil_trt_spp=='uninvaded_uninvaded'|soil_trt_spp=='untreated_GEMO')))

#cysc invasion
MPNcysc <- MPNinv%>%
  select(species, soil_site, uninvaded_uninvaded, untreated_CYSC)%>%
  filter(species!='SPJU', species!='GEMO')%>%
  gather(key=soil_trt_spp, value=ln_MPN, uninvaded_uninvaded:untreated_CYSC, na.rm=T)

cyscModel <- lmer(ln_MPN ~ species*soil_trt_spp + (1 | soil_site), data=MPNcysc)
summary(cyscModel)
ranef(cyscModel)
fixef(cyscModel)

summary(cyscModelglm <- glm(ln_MPN ~ species*soil_trt_spp, data=MPNcysc))

#spju invasion
MPNspju <- MPNinv%>%
  select(species, soil_site, uninvaded_uninvaded, untreated_SPJU)%>%
  filter(species!='CYSC', species!='GEMO')%>%
  gather(key=soil_trt_spp, value=ln_MPN, uninvaded_uninvaded:untreated_SPJU, na.rm=T)

spjuModel <- lmer(ln_MPN ~ species*soil_trt_spp + (1 | soil_site), data=MPNspju)
summary(spjuModel)
ranef(spjuModel)
fixef(spjuModel)

summary(spjuModelglm <- glm(ln_MPN ~ species*soil_trt_spp, data=MPNspju))

#gemo removal
MPNgemo <- MPNinv%>%
  select(species, soil_site, herbicided_GEMO, mowed_GEMO, pulled_GEMO, uninvaded_uninvaded, untreated_GEMO)%>%
  filter(species!='SPJU', species!='CYSC')%>%
  gather(key=soil_trt_spp, value=ln_MPN, herbicided_GEMO:untreated_GEMO, na.rm=T)

gemoRemModel <- lmer(ln_MPN ~ species*soil_trt_spp + (1 | soil_site), data=MPNgemo)
summary(gemoRemModel)
ranef(gemoRemModel)
fixef(gemoRemModel)

summary(gemoRemModelglm <- glm(ln_MPN ~ species*soil_trt_spp, data=MPNgemo))


###figures of invaded vs uninvaded areas

#bar graphs - nice for showing pattern of significance
cyscEffectPlot <- ggplot(data=barGraphStats(data=MPNcysc, variable="ln_MPN", byFactorNames=c("species", "soil_trt_spp")), aes(x=species, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('') +
  theme(legend.position=c(0.9,0.8)) +
  scale_y_continuous(limits=c(0,8)) +
  scale_x_discrete(limits=c('CYSC', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

gemoEffectPlot <- ggplot(data=barGraphStats(data=subset(MPNgemo, soil_trt_spp=='uninvaded_uninvaded'|soil_trt_spp=='untreated_GEMO'), variable="ln_MPN", byFactorNames=c("species", "soil_trt_spp")), aes(x=species, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('ln Most Probable Number (MPN)') +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,8)) +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

spjuEffectPlot <- ggplot(data=barGraphStats(data=MPNspju, variable="ln_MPN", byFactorNames=c("species", "soil_trt_spp")), aes(x=species, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('') +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,8)) +
  scale_x_discrete(limits=c('SPJU', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

#3 panel figure
pushViewport(viewport(layout=grid.layout(3,1)))
print(cyscEffectPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(spjuEffectPlot, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(gemoEffectPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))


#boxplots - nice for showing that the CYSC pattern is mainly driven by one high point at blithedale summit
cyscEffectPlot <- ggplot(data=MPNcysc, aes(x=species, y=ln_MPN, fill=soil_trt_spp)) +
  geom_boxplot() + 
  xlab('') +
  ylab('') +
  theme(legend.position=c(0.1,0.85)) +
  scale_y_continuous(limits=c(0,8)) +
  scale_x_discrete(limits=c('CYSC', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

gemoEffectPlot <- ggplot(data=subset(MPNgemo, soil_trt_spp=='uninvaded_uninvaded'|soil_trt_spp=='untreated_GEMO'), aes(x=species, y=ln_MPN, fill=soil_trt_spp)) +
  geom_boxplot() +
  xlab('') +
  ylab('ln MPN per g soil') +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,8)) +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

spjuEffectPlot <- ggplot(data=MPNspju, aes(x=species, y=ln_MPN, fill=soil_trt_spp)) +
  geom_boxplot() +
  xlab('Plant Species') +
  ylab('') +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,8)) +
  scale_x_discrete(limits=c('SPJU', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

#3 panel figure
pushViewport(viewport(layout=grid.layout(3,1)))
print(cyscEffectPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(spjuEffectPlot, vp=viewport(layout.pos.row=3, layout.pos.col=1)) 
print(gemoEffectPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))



#species responses to GEMO invasions and treatments
#get proper treatment order
MPNgemoOrder <- MPNgemo%>%
  mutate(order=ifelse(soil_trt_spp=='uninvaded_uninvaded', 'a', ifelse(soil_trt_spp=='untreated_GEMO', 'b', ifelse(soil_trt_spp=='pulled_GEMO', 'c', ifelse(soil_trt_spp=='herbicided_GEMO', 'd', 'e')))))

#bar graph
ggplot(data=barGraphStats(data=MPNgemoOrder, variable="ln_MPN", byFactorNames=c("species", "order")), aes(x=species, y=mean, fill=order)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('ln MPN per g soil') +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  theme(legend.position=c(0.9,0.9), axis.title.y=element_text(margin=margin(r=10))) +
  scale_fill_manual(labels=c('uninvaded', 'invaded', 'pulled', 'herbicide', 'mowed'), values=c('#009900', '#FF9900', '#A8A9FF', '#7A75CE', '#554C9E'))

#boxplot
ggplot(data=arrange(MPNgemoOrder, order), aes(x=species, y=ln_MPN, fill=order)) +
  geom_boxplot() +
  xlab('Plant Species') +
  ylab('ln MPN per g soil') +
  # theme(legend.position='none') +
  # scale_y_continuous(limits=c(0,8)) +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded', 'pulled', 'herbicide', 'mowed'), values=c('#009900', '#FF9900', '#A8A9FF', '#7A75CE', '#554C9E'))








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

