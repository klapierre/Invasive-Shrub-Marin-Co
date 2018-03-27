library(ggplot2)
library(grid)
library(tidyverse)
library(lme4)
library(lsmeans)

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
#invasion models
acglModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='ACGL'))
summary(acglModel)

lunaModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='LUNA'))
summary(lunaModel)

cyscModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='CYSC'))
summary(cyscModel)

gemoModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='GEMO'))
summary(gemoModel)

spjuModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='SPJU'))
summary(spjuModel)


###figures of invaded vs uninvaded areas

#bar graphs - nice for showing pattern of significance
acglEffectPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='ACGL'), variable="total_mass", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('') +
  theme(legend.position='none') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_CYSC', 'untreated_GEMO', 'untreated_SPJU'), labels=c('uninvaded', 'CYSC', 'GEMO', 'SPJU')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded', 'invaded', 'invaded'), values=c("#009900", "#FF9900", "#FF9900", "#FF9900"))  +
  annotate('text', x=1, y=0.024, label='a', size=12) +
  annotate('text', x=2, y=0.036, label='b', size=12) +
  annotate('text', x=3, y=0.031, label='b', size=12) +
  annotate('text', x=4, y=0.032, label='ab', size=12) +
  annotate('text', x=0.5, y=0.04, label='(a) ACGL', size=12, hjust='left')

lunaEffectPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='LUNA'), variable="total_mass", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('') +
  theme(legend.position='none') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_CYSC', 'untreated_GEMO', 'untreated_SPJU'), labels=c('uninvaded', 'CYSC', 'GEMO', 'SPJU')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded', 'invaded', 'invaded'), values=c("#009900", "#FF9900", "#FF9900", "#FF9900"))  +
  annotate('text', x=1, y=0.53, label='a', size=12) +
  annotate('text', x=2, y=0.92, label='b', size=12) +
  annotate('text', x=3, y=0.65, label='ab', size=12) +
  annotate('text', x=4, y=0.83, label='b', size=12) +
  annotate('text', x=0.5, y=1, label='(b) LUNA', size=12, hjust='left')

cyscEffectPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='CYSC'), variable="total_mass", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('') +
  theme(legend.position='none') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_CYSC', 'untreated_GEMO', 'untreated_SPJU'), labels=c('uninvaded', 'CYSC', 'GEMO', 'SPJU')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded', 'invaded', 'invaded'), values=c("#009900", "#FF9900", "#FF9900", "#FF9900")) +
  annotate('text', x=0.5, y=0.13, label='(c) CYSC', size=12, hjust='left')

gemoEffectPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='GEMO'), variable="total_mass", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('') +
  theme(legend.position='none') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_CYSC', 'untreated_GEMO', 'untreated_SPJU'), labels=c('uninvaded', 'CYSC', 'GEMO', 'SPJU')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded', 'invaded', 'invaded'), values=c("#009900", "#FF9900", "#FF9900", "#FF9900"))  +
  annotate('text', x=1, y=0.081, label='a', size=12) +
  annotate('text', x=2, y=0.11, label='b', size=12) +
  annotate('text', x=3, y=0.13, label='bc', size=12) +
  annotate('text', x=4, y=0.153, label='c', size=12) +
  annotate('text', x=0.5, y=0.153, label='(d) GEMO', size=12, hjust='left')

spjuEffectPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp!='NA_NA'&soil_trt_spp!='mowed_GEMO'&soil_trt_spp!='pulled_GEMO'&soil_trt_spp!='herbicided_GEMO'&species=='SPJU'), variable="total_mass", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  xlab('') +
  ylab('') +
  theme(legend.position='none') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_CYSC', 'untreated_GEMO', 'untreated_SPJU'), labels=c('uninvaded', 'CYSC', 'GEMO', 'SPJU')) +
  scale_fill_manual(labels=c('uninvaded', 'invaded', 'invaded', 'invaded'), values=c("#009900", "#FF9900", "#FF9900", "#FF9900"))  +
  annotate('text', x=1, y=0.165, label='a', size=12) +
  annotate('text', x=2, y=0.195, label='b', size=12) +
  annotate('text', x=3, y=0.194, label='b', size=12) +
  annotate('text', x=4, y=0.21, label='b', size=12) +
  annotate('text', x=0.5, y=0.21, label='(e) SPJU', size=12, hjust='left')

#5 panel figure
pushViewport(viewport(layout=grid.layout(2,3)))
print(acglEffectPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lunaEffectPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(cyscEffectPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(gemoEffectPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(spjuEffectPlot, vp=viewport(layout.pos.row=2, layout.pos.col=3))
#export at 1500x2000


# #proportion difference
# uninvaded <- biomass%>%
#   filter(soil_trt=='uninvaded')%>%
#   group_by(species, soil_site)%>%
#   summarise(total_mass_uninv=mean(total_mass), root_uninv=mean(root), shoot_uninv=mean(shoot))%>%
#   ungroup()
# 
#   
# proportionDifference <- biomass%>%
#   filter(soil_trt!='uninvaded')%>%
#   left_join(uninvaded)%>%
#   mutate(total_diff=(total_mass-total_mass_uninv)/total_mass_uninv, root_diff=(root-root_uninv)/root_uninv, shoot_diff=(shoot-shoot_uninv)/shoot_uninv)
# 
# 
# #gemo removal
# acglRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='ACGL'))
# summary(acglRemModel)
# lsmeans(acglRemModel, pairwise~soil_trt_spp)
# 
# lunaRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='LUNA'))
# summary(lunaRemModel)
# lsmeans(lunaRemModel, pairwise~soil_trt_spp)
# 
# cyscRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='CYSC'))
# summary(cyscRemModel)
# lsmeans(cyscRemModel, pairwise~soil_trt_spp)
# 
# gemoRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='GEMO'))
# summary(gemoRemModel)
# lsmeans(gemoRemModel, pairwise~soil_trt_spp)
# 
# spjuRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='SPJU'))
# summary(spjuRemModel)
# lsmeans(spjuRemModel, pairwise~soil_trt_spp)
# 
# #species responses to GEMO invasions and treatments
# #bar graph
# acglRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='ACGL'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
#   ylab('Proportional Biomass Difference') +
#   scale_x_discrete(limits=c('untreated_GEMO', 'pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('untrt.', 'pulled', 'herbic.', 'mowed')) +
#   theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
#         axis.text.x  = element_text(angle=90, vjust=0.5))  +
#   scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF', '#FF9900')) +
#   # annotate('text', x=1, y=0.165, label='a', size=12) +
#   # annotate('text', x=2, y=0.195, label='b', size=12) +
#   # annotate('text', x=3, y=0.194, label='b', size=12) +
#   # annotate('text', x=4, y=0.21, label='b', size=12) +
#   annotate('text', x=0.5, y=1.9, label='(a) ACGL', size=8, hjust='left')
# 
# lunaRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='LUNA'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
#   ylab('Proportional Biomass Difference') +
#   scale_x_discrete(limits=c('untreated_GEMO', 'pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('untrt.', 'pulled', 'herbic.', 'mowed')) +
#   theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
#         axis.text.x  = element_text(angle=90, vjust=0.5))  +
#   scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF', '#FF9900')) +
#   annotate('text', x=1, y=1.9, label='ab', size=8) +
#   annotate('text', x=2, y=1.45, label='a', size=8) +
#   annotate('text', x=3, y=4.3, label='b', size=8) +
#   annotate('text', x=4, y=2.75, label='b', size=8) +
#   annotate('text', x=0.5, y=4.2, label='(b) LUNA', size=8, hjust='left')
# 
# cyscRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='CYSC'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
#   ylab('Proportional Biomass Difference') +
#   scale_x_discrete(limits=c('untreated_GEMO', 'pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('untrt.', 'pulled', 'herbic.', 'mowed')) +
#   theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
#         axis.text.x  = element_text(angle=90, vjust=0.5))  +
#   scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF', '#FF9900')) +
#   # annotate('text', x=1, y=0.8, label='a', size=8) +
#   # annotate('text', x=2, y=1.3, label='ab', size=8) +
#   # annotate('text', x=3, y=1.8, label='b', size=8) +
#   # annotate('text', x=4, y=0.7, label='a', size=8) +
#   annotate('text', x=0.5, y=1.7, label='(c) CYSC', size=8, hjust='left')
# 
# gemoRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='GEMO'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
#   ylab('Proportional Biomass Difference') +
#   scale_x_discrete(limits=c('untreated_GEMO', 'pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('untrt.', 'pulled', 'herbic.', 'mowed')) +
#   theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
#         axis.text.x  = element_text(angle=90, vjust=0.5))  +
#   scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF', '#FF9900')) +
#   annotate('text', x=1, y=1.5, label='a', size=8) +
#   annotate('text', x=2, y=0.9, label='b', size=8) +
#   annotate('text', x=3, y=1.2, label='ab', size=8) +
#   annotate('text', x=4, y=1.1, label='ab', size=8) +
#   annotate('text', x=0.5, y=1.7, label='(d) GEMO', size=8, hjust='left')
# 
# spjuRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp!='NA_NA'&soil_trt_spp!='untreated_SPJU'&soil_trt_spp!='untreated_CYSC'&species=='SPJU'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
#   ylab('Proportional Biomass Difference') +
#   scale_x_discrete(limits=c('untreated_GEMO', 'pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('untrt.', 'pulled', 'herbic.', 'mowed')) +
#   theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
#         axis.text.x  = element_text(angle=90, vjust=0.5))  +
#   scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF', '#FF9900')) +
#   # annotate('text', x=1, y=0.165, label='a', size=12) +
#   # annotate('text', x=2, y=0.195, label='b', size=12) +
#   # annotate('text', x=3, y=0.194, label='b', size=12) +
#   # annotate('text', x=4, y=0.21, label='b', size=12) +
#   annotate('text', x=0.5, y=0.6, label='(e) SPJU', size=8, hjust='left')
# 
# #5 panel figure
# pushViewport(viewport(layout=grid.layout(2,3)))
# print(acglRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(lunaRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(cyscRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
# print(gemoRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
# print(spjuRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=3))
# #export at 1500x2000



###compare uninvaded vs invaded

#gemo removal
acglRemModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='ACGL'))
summary(acglRemModel)
lsmeans(acglRemModel, pairwise~soil_trt_spp)

lunaRemModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='LUNA'))
summary(lunaRemModel)
lsmeans(lunaRemModel, pairwise~soil_trt_spp)

cyscRemModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='CYSC'))
summary(cyscRemModel)
lsmeans(cyscRemModel, pairwise~soil_trt_spp)

gemoRemModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='GEMO'))
summary(gemoRemModel)
lsmeans(gemoRemModel, pairwise~soil_trt_spp)

spjuRemModel <- lmer(total_mass ~ soil_trt_spp + (1 | soil_site), data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='SPJU'))
summary(spjuRemModel)
lsmeans(spjuRemModel, pairwise~soil_trt_spp)

#species responses to GEMO invasions and treatments
#bar graph
acglRemovalPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='ACGL'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Total Biomass (g)') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#009900', '#FF9900')) +
  annotate('text', x=1, y=0.026, label='a', size=8) +
  annotate('text', x=2, y=0.032, label='b', size=8) +
  annotate('text', x=0.5, y=0.035, label='(a) ACGL', size=8, hjust='left')

lunaRemovalPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='LUNA'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Total Biomass (g)') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#009900', '#FF9900')) +
  annotate('text', x=1, y=0.56, label='a', size=8) +
  annotate('text', x=2, y=0.68, label='b', size=8) +
  annotate('text', x=0.5, y=0.7, label='(b) LUNA', size=8, hjust='left')

cyscRemovalPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='CYSC'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Total Biomass (g)') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#009900', '#FF9900')) +
  # annotate('text', x=1, y=0.8, label='a', size=8) +
  # annotate('text', x=2, y=1.3, label='ab', size=8) +
  annotate('text', x=0.5, y=0.13, label='(c) CYSC', size=8, hjust='left')

gemoRemovalPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='GEMO'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Total Biomass (g)') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#009900', '#FF9900')) +
  annotate('text', x=1, y=0.09, label='a', size=8) +
  annotate('text', x=2, y=0.14, label='b', size=8) +
  annotate('text', x=0.5, y=0.15, label='(d) GEMO', size=8, hjust='left')

spjuRemovalPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='SPJU'), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Total Biomass (g)') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#009900', '#FF9900')) +
  annotate('text', x=1, y=0.17, label='a', size=8) +
  annotate('text', x=2, y=0.2, label='b', size=8) +
  annotate('text', x=0.5, y=0.21, label='(e) SPJU', size=8, hjust='left')

#5 panel figure
pushViewport(viewport(layout=grid.layout(2,3)))
print(acglRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lunaRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(cyscRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(gemoRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(spjuRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=3))
#export at 1500x2000



###proportion difference between removal and invaded (untreated)
uninvaded <- biomass%>%
  filter(soil_trt=='uninvaded')%>%
  group_by(species, soil_site)%>%
  summarise(total_mass_uninv=mean(total_mass), root_uninv=mean(root), shoot_uninv=mean(shoot))%>%
  ungroup()


proportionDifference <- biomass%>%
  filter(soil_trt!='untreated_GEMO')%>%
  left_join(uninvaded)%>%
  mutate(total_diff=(total_mass-total_mass_uninv)/total_mass_uninv, root_diff=(root-root_uninv)/root_uninv, shoot_diff=(shoot-shoot_uninv)/shoot_uninv)


#gemo removal
acglRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='ACGL'))
summary(acglRemModel)
lsmeans(acglRemModel, pairwise~soil_trt_spp)
#t-tests to identify differences from 1.0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(proportionDifference, soil_trt_spp=='pulled_GEMO'&species=='ACGL')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='herbicided_GEMO'&species=='ACGL')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='mowed_GEMO'&species=='ACGL')
t.test(test$total_diff, mu=1)

lunaRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='LUNA'))
summary(lunaRemModel)
lsmeans(lunaRemModel, pairwise~soil_trt_spp)
#t-tests to identify differences from 1.0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(proportionDifference, soil_trt_spp=='pulled_GEMO'&species=='LUNA')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='herbicided_GEMO'&species=='LUNA')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='mowed_GEMO'&species=='LUNA')
t.test(test$total_diff, mu=1)

cyscRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='CYSC'))
summary(cyscRemModel)
lsmeans(cyscRemModel, pairwise~soil_trt_spp)
#t-tests to identify differences from 1.0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(proportionDifference, soil_trt_spp=='pulled_GEMO'&species=='CYSC')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='herbicided_GEMO'&species=='CYSC')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='mowed_GEMO'&species=='CYSC')
t.test(test$total_diff, mu=1)

gemoRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='GEMO'))
summary(gemoRemModel)
lsmeans(gemoRemModel, pairwise~soil_trt_spp)
#t-tests to identify differences from 1.0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(proportionDifference, soil_trt_spp=='pulled_GEMO'&species=='GEMO')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='herbicided_GEMO'&species=='GEMO')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='mowed_GEMO'&species=='GEMO')
t.test(test$total_diff, mu=1)

spjuRemModel <- lmer(total_diff ~ soil_trt_spp + (1 | soil_site), data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='SPJU'))
summary(spjuRemModel)
lsmeans(spjuRemModel, pairwise~soil_trt_spp)
#t-tests to identify differences from 1.0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(proportionDifference, soil_trt_spp=='pulled_GEMO'&species=='SPJU')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='herbicided_GEMO'&species=='SPJU')
t.test(test$total_diff, mu=1)
test <- subset(proportionDifference, soil_trt_spp=='mowed_GEMO'&species=='SPJU')
t.test(test$total_diff, mu=1)

#species responses to GEMO invasions and treatments
#bar graph
acglRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='ACGL'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Proportional Biomass Difference') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF')) +
  geom_abline(intercept=1, slope=0, linetype=2) +
  annotate('text', x=2, y=0.6, label='*', size=8) +
  annotate('text', x=0.5, y=1.9, label='(a) ACGL', size=8, hjust='left')

lunaRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='LUNA'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Proportional Biomass Difference') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF')) +
  geom_abline(intercept=1, slope=0, linetype=2) +
  annotate('text', x=1, y=1.4, label='a', size=8) +
  annotate('text', x=2, y=4.2, label='b', size=8) +
  annotate('text', x=3, y=2.7, label='b', size=8) +
  annotate('text', x=0.5, y=4.2, label='(b) LUNA', size=8, hjust='left')

cyscRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='CYSC'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Proportional Biomass Difference') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF')) +
  geom_abline(intercept=1, slope=0, linetype=2) +
  annotate('text', x=1, y=1.3, label='ab', size=8) +
  annotate('text', x=2, y=1.8, label='a', size=8) +
  annotate('text', x=3, y=0.7, label='b*', size=8) +
  annotate('text', x=0.5, y=1.8, label='(c) CYSC', size=8, hjust='left')

gemoRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='GEMO'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Proportional Biomass Difference') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF')) +
  geom_abline(intercept=1, slope=0, linetype=2) +
  annotate('text', x=1, y=0.77, label='*', size=8) +
  annotate('text', x=0.5, y=1.2, label='(d) GEMO', size=8, hjust='left')

spjuRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='SPJU'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Proportional Biomass Difference') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5))  +
  scale_fill_manual(values=c('#7A75CE', '#554C9E', '#A8A9FF')) +
  geom_abline(intercept=1, slope=0, linetype=2) +
  annotate('text', x=1, y=0.3, label='*', size=8) +
  annotate('text', x=2, y=0.31, label='*', size=8) +
  annotate('text', x=3, y=0.61, label='*', size=8) +
  annotate('text', x=0.5, y=1.1, label='(e) SPJU', size=8, hjust='left')

#5 panel figure
pushViewport(viewport(layout=grid.layout(2,3)))
print(acglRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lunaRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(cyscRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(gemoRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(spjuRemovalPlot, vp=viewport(layout.pos.row=2, layout.pos.col=3))
#export at 1500x2000