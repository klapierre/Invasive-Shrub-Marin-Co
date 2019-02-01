library(lme4)
library(lmerTest)
library(grid)
library(tidyverse)

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


#MPN data
trt <- read.csv('La Pierre_MPN_treatments_2015.csv')
pots <- read.csv('La Pierre_MPN_pots_2015.csv')
nods <- read.csv('La Pierre_MPN_harvest_2015.csv')%>%
  filter(notes!='no growth', notes!='dead')%>%
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
#comparing MPN invaded vs uninvaded
MPNgemo <- MPNinv%>%
  select(species, soil_site, uninvaded_uninvaded, untreated_GEMO)%>%
  filter(untreated_GEMO!='NA')%>%
  gather(key=soil_trt_spp, value=ln_MPN, uninvaded_uninvaded:untreated_GEMO, na.rm=T)

summary(acglMPNModel <- lmer(ln_MPN ~ soil_trt_spp + (1|soil_site), data=subset(MPNgemo, species=='ACGL')))
summary(acglMPNModel <- lmer(ln_MPN ~ soil_trt_spp + (1|soil_site), data=subset(MPNgemo, species=='LUBI')))
summary(acglMPNModel <- lmer(ln_MPN ~ soil_trt_spp + (1|soil_site), data=subset(MPNgemo, species=='LUNA')))
summary(acglMPNModel <- lmer(ln_MPN ~ soil_trt_spp + (1|soil_site), data=subset(MPNgemo, species=='GEMO')))

#bar graph
acglRemovalPlot <- ggplot(data=barGraphStats(data=subset(MPNgemo, species=='ACGL'), variable="ln_MPN", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('ln MPN') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080')) +
  annotate('text', x=1, y=0.15, label='a', size=8) +
  annotate('text', x=2, y=2.5, label='b', size=8) +
  annotate('text', x=0.5, y=2.5, label=expression(paste('(a) ',italic('A. glaber'))), size=8, hjust='left')

lubiRemovalPlot <- ggplot(data=barGraphStats(data=subset(MPNgemo, species=='LUBI'), variable="ln_MPN", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080')) +
  # annotate('text', x=1, y=0.15, label='a', size=8) +
  # annotate('text', x=2, y=3.2, label='b', size=8) +
  annotate('text', x=0.5, y=3.2, label=expression(paste('(b) ',italic('L. bicolor'))), size=8, hjust='left')

lunaRemovalPlot <- ggplot(data=barGraphStats(data=subset(MPNgemo, species=='LUNA'), variable="ln_MPN", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080')) +
  annotate('text', x=1, y=0.6, label='a', size=8) +
  annotate('text', x=2, y=4, label='b', size=8) +
  annotate('text', x=0.5, y=4, label=expression(paste('(c) ',italic('L. nanus'))), size=8, hjust='left')

gemoRemovalPlot <- ggplot(data=barGraphStats(data=subset(MPNgemo, species=='GEMO'), variable="ln_MPN", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080')) +
  annotate('text', x=1, y=0.7, label='a', size=8) +
  annotate('text', x=2, y=5.5, label='b', size=8) +
  annotate('text', x=0.5, y=6.0, label=expression(paste('(d) ',italic('G. monspessulana'))), size=8, hjust='left')

#4 panel figure
pushViewport(viewport(layout=grid.layout(1,4)))
print(acglRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lubiRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(lunaRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(gemoRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=4))
#export at 2000x600



#comparing MPN with gemo removal techniques (proportional to invaded MPN)
MPNuntreated <- MPNinv%>%
  select(species, soil_site, untreated_GEMO)%>%
  filter(untreated_GEMO!='NA')

MPNgemoRem <- MPNinv%>%
  select(species, soil_site, herbicided_GEMO, mowed_GEMO, pulled_GEMO)%>%
  gather(key=soil_trt_spp, value=ln_MPN, herbicided_GEMO:pulled_GEMO, na.rm=T)%>%
  left_join(MPNuntreated)%>%
  mutate(proportion_MPN=((ln_MPN+0.1)-(untreated_GEMO+0.1))/(untreated_GEMO+0.1))

summary(aov(proportion_MPN ~ soil_trt_spp, data=subset(MPNgemoRem, species=='ACGL')))
#t-tests to identify differences from 0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(MPNgemoRem, soil_trt_spp=='pulled_GEMO'&species=='ACGL')
t.test(test$proportion_MPN, mu=0)
test <- subset(MPNgemoRem, soil_trt_spp=='herbicided_GEMO'&species=='ACGL')
t.test(test$proportion_MPN, mu=0)
test <- subset(MPNgemoRem, soil_trt_spp=='mowed_GEMO'&species=='ACGL')
t.test(test$proportion_MPN, mu=0)

summary(aov(proportion_MPN ~ soil_trt_spp, data=subset(MPNgemoRem, species=='LUBI')))
#t-tests to identify differences from 0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(MPNgemoRem, soil_trt_spp=='pulled_GEMO'&species=='LUBI')
t.test(test$proportion_MPN, mu=0)
test <- subset(MPNgemoRem, soil_trt_spp=='herbicided_GEMO'&species=='LUBI')
t.test(test$proportion_MPN, mu=0)
test <- subset(MPNgemoRem, soil_trt_spp=='mowed_GEMO'&species=='LUBI')
t.test(test$proportion_MPN, mu=0)

summary(aov(proportion_MPN ~ soil_trt_spp, data=subset(MPNgemoRem, species=='LUNA')))
#t-tests to identify differences from 0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(MPNgemoRem, soil_trt_spp=='pulled_GEMO'&species=='LUNA')
t.test(test$proportion_MPN, mu=0)
test <- subset(MPNgemoRem, soil_trt_spp=='herbicided_GEMO'&species=='LUNA')
t.test(test$proportion_MPN, mu=0)
test <- subset(MPNgemoRem, soil_trt_spp=='mowed_GEMO'&species=='LUNA')
t.test(test$proportion_MPN, mu=0)

summary(aov(proportion_MPN ~ soil_trt_spp, data=subset(MPNgemoRem, species=='GEMO')))
#t-tests to identify differences from 0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(MPNgemoRem, soil_trt_spp=='pulled_GEMO'&species=='GEMO')
t.test(test$proportion_MPN, mu=0)
test <- subset(MPNgemoRem, soil_trt_spp=='herbicided_GEMO'&species=='GEMO')
t.test(test$proportion_MPN, mu=0)
test <- subset(MPNgemoRem, soil_trt_spp=='mowed_GEMO'&species=='GEMO')
t.test(test$proportion_MPN, mu=0)


#bar graph
acglRemovalPlot <- ggplot(data=barGraphStats(data=subset(MPNgemoRem, species=='ACGL'), variable="proportion_MPN", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Proportional MPN Difference') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank()) +
  geom_abline(intercept=0, slope=0, linetype=1) +
  # annotate('text', x=1, y=-1.5, label='*', size=8) +
  annotate('text', x=0.5, y=6, label='(a) ACGL', size=8, hjust='left')

lubiRemovalPlot <- ggplot(data=barGraphStats(data=subset(MPNgemoRem, species=='LUBI'), variable="proportion_MPN", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank()) +
  geom_abline(intercept=0, slope=0, linetype=1) +
  # annotate('text', x=1, y=-3, label='*', size=8) +
  annotate('text', x=0.5, y=30, label='(b) LUBI', size=8, hjust='left')

lunaRemovalPlot <- ggplot(data=barGraphStats(data=subset(MPNgemoRem, species=='LUNA'), variable="proportion_MPN", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank()) +
  geom_abline(intercept=0, slope=0, linetype=1) +
  # annotate('text', x=1, y=-1.1, label='*', size=8) +
  # annotate('text', x=2, y=-0.75, label='*', size=8) +
  # annotate('text', x=3, y=0.4, label='*', size=8) +
  annotate('text', x=0.5, y=1.5, label='(c) LUNA', size=8, hjust='left')

gemoRemovalPlot <-  ggplot(data=barGraphStats(data=subset(MPNgemoRem, species=='GEMO'), variable="proportion_MPN", byFactorNames=c("soil_trt_spp")), aes(x=soil_trt_spp, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank()) +
  geom_abline(intercept=0, slope=0, linetype=1) +
  annotate('text', x=1, y=-1.1, label='*', size=8) +
  # annotate('text', x=2, y=-0.75, label='*', size=8) +
  # annotate('text', x=3, y=0.4, label='*', size=8) +
  annotate('text', x=0.5, y=1.5, label='(d) GEMO', size=8, hjust='left')

#4 panel figure
pushViewport(viewport(layout=grid.layout(1,4)))
print(acglRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lubiRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(lunaRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(gemoRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=4))
#export at 2000x600





###biomass - whole soil inoculations

#data
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

###compare uninvaded vs invaded
#gemo removal
summary(acglRemModel <- aov(total_mass ~ soil_trt_spp, data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='ACGL'&!(soil_site %in% c('boyd memorial', 'deer park', 'horse hill', 'larkspur ferry', 'loma alta', 'phoenix lake')))))
summary(lunaRemModel <- aov(total_mass ~ soil_trt_spp, data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='LUNA'&!(soil_site %in% c('boyd memorial', 'deer park', 'horse hill', 'larkspur ferry', 'loma alta', 'phoenix lake')))))
summary(gemoRemModel <- aov(total_mass ~ soil_trt_spp, data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='GEMO'&!(soil_site %in% c('boyd memorial', 'deer park', 'horse hill', 'larkspur ferry', 'loma alta', 'phoenix lake')))))

#species responses to GEMO invasions and treatments
#bar graph
acglRemovalPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='ACGL'&!(soil_site %in% c('boyd memorial', 'deer park', 'horse hill', 'larkspur ferry', 'loma alta', 'phoenix lake'))), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Total Biomass (g)') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080')) +
  annotate('text', x=1, y=0.02, label='a', size=8) +
  annotate('text', x=2, y=0.032, label='b', size=8) +
  annotate('text', x=0.5, y=0.032, label='(a) ACGL', size=8, hjust='left')

lunaRemovalPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='LUNA'&!(soil_site %in% c('boyd memorial', 'deer park', 'horse hill', 'larkspur ferry', 'loma alta', 'phoenix lake'))), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080')) +
  annotate('text', x=1, y=0.39, label='a', size=8) +
  annotate('text', x=2, y=0.67, label='b', size=8) +
  annotate('text', x=0.5, y=0.67, label='(b) LUNA', size=8, hjust='left')

gemoRemovalPlot <- ggplot(data=barGraphStats(data=subset(biomass, soil_trt_spp %in% c('untreated_GEMO', 'uninvaded_uninvaded')&species=='GEMO'&!(soil_site %in% c('boyd memorial', 'deer park', 'horse hill', 'larkspur ferry', 'loma alta', 'phoenix lake'))), variable="total_mass", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('uninvaded_uninvaded', 'untreated_GEMO'), labels=c('uninv.', 'untrt.')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080')) +
  annotate('text', x=1, y=0.07, label='a', size=8) +
  annotate('text', x=2, y=0.135, label='b', size=8) +
  annotate('text', x=0.5, y=0.14, label='(c) GEMO', size=8, hjust='left')

#3 panel figure
pushViewport(viewport(layout=grid.layout(1,3)))
print(acglRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lunaRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(gemoRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export at 1500x600



###proportion difference between removal and invaded (untreated)
untreated <- biomass%>%
  filter(soil_trt=='untreated', !(soil_site %in% c('boyd memorial', 'deer park', 'horse hill', 'larkspur ferry', 'loma alta', 'phoenix lake')))%>%
  group_by(species, soil_site)%>%
  summarise(total_mass_untrt=mean(total_mass), root_untrt=mean(root), shoot_untrt=mean(shoot))%>%
  ungroup()

proportionDifference <- biomass%>%
  filter(soil_trt!='untreated', soil_trt!='uninvaded', !(soil_site %in% c('boyd memorial', 'deer park', 'horse hill', 'larkspur ferry', 'loma alta', 'phoenix lake')))%>%
  left_join(untreated)%>%
  mutate(total_diff=(total_mass-total_mass_untrt)/total_mass_untrt, root_diff=(root-root_untrt)/root_untrt, shoot_diff=(shoot-shoot_untrt)/shoot_untrt)

#gemo removal
summary(acglRemModel <- aov(total_diff ~ soil_trt_spp, data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='ACGL')))
#t-tests to identify differences from 0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(proportionDifference, soil_trt_spp=='pulled_GEMO'&species=='ACGL')
t.test(test$total_diff, mu=0)
test <- subset(proportionDifference, soil_trt_spp=='herbicided_GEMO'&species=='ACGL')
t.test(test$total_diff, mu=0)
test <- subset(proportionDifference, soil_trt_spp=='mowed_GEMO'&species=='ACGL')
t.test(test$total_diff, mu=0)

summary(lunaRemModel <- aov(total_diff ~ soil_trt_spp, data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='LUNA')))
#t-tests to identify differences from 0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(proportionDifference, soil_trt_spp=='pulled_GEMO'&species=='LUNA')
t.test(test$total_diff, mu=0)
test <- subset(proportionDifference, soil_trt_spp=='herbicided_GEMO'&species=='LUNA')
t.test(test$total_diff, mu=0)
test <- subset(proportionDifference, soil_trt_spp=='mowed_GEMO'&species=='LUNA')
t.test(test$total_diff, mu=0)

summary(gemoRemModel <- aov(total_diff ~ soil_trt_spp, data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='GEMO')))
#t-tests to identify differences from 0 (no proportional difference from untreated GEMO invaded areas)
test <- subset(proportionDifference, soil_trt_spp=='pulled_GEMO'&species=='GEMO')
t.test(test$total_diff, mu=0)
test <- subset(proportionDifference, soil_trt_spp=='herbicided_GEMO'&species=='GEMO')
t.test(test$total_diff, mu=0)
test <- subset(proportionDifference, soil_trt_spp=='mowed_GEMO'&species=='GEMO')
t.test(test$total_diff, mu=0)

#species responses to GEMO invasions and treatments
#bar graph
acglRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='ACGL'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Proportional Biomass Difference') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080', '#808080')) +
  geom_abline(intercept=0, slope=0, linetype=1) +
  # annotate('text', x=2, y=0.6, label='*', size=8) +
  annotate('text', x=0.5, y=0.3, label='(a) ACGL', size=8, hjust='left')

lunaRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='LUNA'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080', '#808080')) +
  geom_abline(intercept=0, slope=0, linetype=1) +
  annotate('text', x=1, y=-0.7, label='a', size=8) +
  annotate('text', x=2, y=0.82, label='b', size=8) +
  annotate('text', x=3, y=0.75, label='b', size=8) +
  annotate('text', x=0.5, y=0.82, label='(b) LUNA', size=8, hjust='left')

gemoRemovalPlot <- ggplot(data=barGraphStats(data=subset(proportionDifference, soil_trt_spp %in% c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO')&species=='GEMO'), variable="total_diff", byFactorNames=c("species", "soil_trt_spp")), aes(x=soil_trt_spp, y=mean, fill=soil_trt_spp)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=position_dodge(0.9)) +
  ylab('Proportional Biomass Difference') +
  scale_x_discrete(limits=c('pulled_GEMO', 'herbicided_GEMO', 'mowed_GEMO'), labels=c('pulled', 'herbic.', 'mowed')) +
  theme(legend.position='none', axis.title.y=element_text(margin=margin(r=10)), axis.title.x=element_blank())  +
  scale_fill_manual(values=c('#808080', '#808080', '#808080')) +
  geom_abline(intercept=0, slope=0, linetype=1) +
  annotate('text', x=1, y=-0.5, label='a*', size=8) +
  annotate('text', x=2, y=-0.25, label='b', size=8) +
  annotate('text', x=3, y=-0.18, label='b', size=8) +
  annotate('text', x=0.5, y=0.2, label='(c) GEMO', size=8, hjust='left')

#3 panel figure
pushViewport(viewport(layout=grid.layout(1,3)))
print(acglRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lunaRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(gemoRemovalPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export at 1500x600
