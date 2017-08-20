library(ggplot2)
library(grid)
library(tidyr)
library(dplyr)
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


setwd('C:\\Users\\Kim\\Dropbox\\Marin County Invasion\\marin parks data')

trt <- read.csv('La Pierre_MPN_treatments_2015.csv')
pots <- read.csv('La Pierre_whole soil_pots_2016.csv')%>%
  select(-species)
harvest <- read.csv('La Pierre_2016_MC_whole soil_harvest.csv')%>%
  # filter(notes!='no growth', notes!='dead')%>%
  select(-f_nods, -n_nods, -nods_collected, -notes, -harvest_date)%>%
  gather(key=root_shoot_env, value=envelope, root_env, shoot_env)%>%
  mutate(root_shoot=ifelse(root_shoot_env=='root_env', 'root', 'shoot'))
biomass <- read.csv('La Pierre_MC_whole soil_biomass_2017.csv')%>%
  left_join(harvest)%>%
  left_join(pots)%>% #note, the plants that died/didn't grow won't be in this dataframe
  left_join(trt)%>%
  mutate(soil_trt_spp=paste(soil_trt, soil_spp, sep='_'))%>%
  #drop 12 individuals for which there is no pot data
  filter(pot!='NA')%>%
  select(pot, mass_g, species, root_shoot, soil_site, soil_spp, soil_trt, soil_trt_spp)%>%
  #calculate total biomass
  spread(key=root_shoot, value=mass_g)%>%
  mutate(total_mass=root+shoot)%>%
  #drop biomass samples without both root and shoot for now (fix this)
  filter(total_mass!='NA')


# #total biomass response ratio -- highest concentration inoculum (dil=0)
# biomassResponse <- biomassTotal%>%
#   filter(dilution!=1&dilution!=2&dilution!=3)%>%
#   group_by(soil_trt_spp, species)%>%
#   summarise(total_mass=mean(total_mass))%>%
#   #compare by trt
#   select(soil_trt_spp, total_mass, species)%>%
#   spread(key=soil_trt_spp, value=total_mass)%>%
#   #relativize to uninvaded
#   mutate(cyscRR=ifelse(untreated_CYSC=='NA', 'NA', (log(untreated_CYSC/uninvaded_uninvaded))))%>%
#   mutate(spjuRR=ifelse(untreated_SPJU=='NA', 'NA', (log(untreated_SPJU/uninvaded_uninvaded))))%>%
#   mutate(gemoRR=ifelse(untreated_GEMO=='NA', 'NA', (log(untreated_GEMO/uninvaded_uninvaded))))%>%
#   mutate(herbRR=ifelse(herbicided_GEMO=='NA', 'NA', (log(herbicided_GEMO/uninvaded_uninvaded))))%>%
#   mutate(mowRR=ifelse(mowed_GEMO=='NA', 'NA', (log(mowed_GEMO/uninvaded_uninvaded))))%>%
#   mutate(pullRR=ifelse(pulled_GEMO=='NA', 'NA', (log(pulled_GEMO/uninvaded_uninvaded))))


###mixed effects models
#gemo invasion
gemoModel <- lmer(total_mass ~ species*soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_CYSC'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species!='CYSC'&species!='SPJU'))
summary(gemoModel)
ranef(gemoModel)
fixef(gemoModel)
# sjp.lmer(gemoModel, type='fe') #doesn't work for some reason

#cysc invasion
cyscModel <- lmer(total_mass ~ species*soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_GEMO'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species!='GEMO'&species!='SPJU'))
summary(cyscModel)
ranef(cyscModel)
fixef(cyscModel)

#spju invasion
spjuModel <- lmer(total_mass ~ species*soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_GEMO'&soil_trt_spp!='untreated_CYSC'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species!='GEMO'&species!='CYSC'))
summary(spjuModel)
ranef(spjuModel)
fixef(spjuModel)

#gemo removal
gemoRemModel <- lmer(total_mass ~ species*soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species!='SPJU'&species!='CYSC'))
summary(gemoRemModel)
ranef(gemoRemModel)
fixef(gemoRemModel)


###figures of invaded vs uninvaded areas

#bar graphs - nice for showing pattern of significance
cyscEffectPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_GEMO'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species!='GEMO'&species!='SPJU'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=species, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('') +
  theme(legend.position=c(0.78,0.92), axis.title.y=element_text(margin=margin(r=10))) +
  scale_y_continuous(limits=c(0,0.25)) +
  scale_x_discrete(limits=c('CYSC', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

gemoEffectPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_CYSC'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species!='CYSC'&species!='SPJU'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=species, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('Total Biomass (g)') +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10))) +
  scale_y_continuous(limits=c(0,0.25)) +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

spjuEffectPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_GEMO'&soil_trt_spp!='untreated_CYSC'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species!='GEMO'&species!='CYSC'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=species, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('') +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10))) +
  scale_y_continuous(limits=c(0,0.25)) +
  scale_x_discrete(limits=c('SPJU', 'ACGL', 'LUBI', 'LUNA')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded'), values=c("#009900", "#FF9900"))

#3 panel figure
pushViewport(viewport(layout=grid.layout(3,1)))
print(cyscEffectPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(spjuEffectPlot, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(gemoEffectPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))


#species responses to GEMO invasions and treatments
#bar graph
ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species!='SPJU'&species!='CYSC'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=species, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('Plant Species') +
  ylab('Total Biomass (g)') +
  scale_x_discrete(limits=c('GEMO', 'ACGL', 'LUBI', 'LUNA')) +
  theme(legend.position=c(0.9,0.9), axis.title.y=element_text(margin=margin(r=10)))