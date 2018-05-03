library(ggplot2)
library(ellipse)
library(grid)
library(tidyverse)

setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\Marin County Invasion\\marin parks data\\DNA sequences\\REP PCR results\\invaded comparison')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))



#read in nodule identifying information
noduleDatabase <- read.delim('all species and treatments_database.txt')

#read in MDS coordinates for each key and merge with database information
MDScoordinates <- read.delim('all species and treatments_MDS coordinates.txt')%>%
  left_join(noduleDatabase)%>%
  select(-X, -Name)



#MDS - each species planted, uninvaded vs GEMO invaded soils
ACGLinvaded <- ggplot(data=subset(MDScoordinates, subset=(plant_species=='ACGL' & soil_spp %in% c('uninvaded', 'GEMO') & soil_trt %in% c('uninvaded', 'untreated'))), 
       (aes(x=MDS1, y=MDS2, color=soil_trt, shape=soil_trt))) +
  geom_point(size=3) +
  stat_ellipse() +
  scale_color_manual(values=c('black', '#A9A9A9')) +
  annotate('text', x=-0.8, y=0.9, label='(a) ACGL', size=8, hjust='left') +
  theme(legend.position=c(0.85, 0.93))

#too few points for a LUBI plot (only 2 nodules from uninvaded areas)

LUNAinvaded <- ggplot(data=subset(MDScoordinates, subset=(plant_species=='LUNA' & soil_spp %in% c('uninvaded', 'GEMO') & soil_trt %in% c('uninvaded', 'untreated'))), 
                      (aes(x=MDS1, y=MDS2, color=soil_trt, shape=soil_trt))) +
  geom_point(size=3) +
  stat_ellipse() +
  scale_color_manual(values=c('black', '#A9A9A9')) +
  annotate('text', x=-0.8, y=0.9, label='(b) LUNA', size=8, hjust='left') +
  theme(legend.position='none')

GEMOinvaded <- ggplot(data=subset(MDScoordinates, subset=(plant_species=='GEMO' & soil_spp %in% c('uninvaded', 'GEMO') & soil_trt %in% c('uninvaded', 'untreated'))), 
                      (aes(x=MDS1, y=MDS2, color=soil_trt, shape=soil_trt))) +
  geom_point(size=3) +
  stat_ellipse() +
  scale_color_manual(values=c('black', '#A9A9A9')) +
  annotate('text', x=-0.8, y=0.9, label='(c) GEMO', size=8, hjust='left') +
  theme(legend.position='none')

#4 panel figure
pushViewport(viewport(layout=grid.layout(1,3)))
print(ACGLinvaded, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(LUNAinvaded, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(GEMOinvaded, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export at 2000x600



#MDS - each species planted, GEMO invaded vs treated soils
ACGLinvaded <- ggplot(data=subset(MDScoordinates, subset=(plant_species=='ACGL' & soil_spp=='GEMO' & soil_trt %in% c('untreated', 'herbicided', 'mowed', 'pulled'))), 
                      (aes(x=MDS1, y=MDS2, shape=soil_trt))) +
  geom_point(size=3) +
  stat_ellipse() +
  annotate('text', x=-0.8, y=0.9, label='(a) ACGL', size=8, hjust='left') +
  theme(legend.position=c(0.85, 0.90))

#too few points for a LUBI plot (only 2 nodules from uninvaded areas)

LUNAinvaded <- ggplot(data=subset(MDScoordinates, subset=(plant_species=='LUNA' & soil_spp=='GEMO' & soil_trt %in% c('untreated', 'herbicided', 'mowed', 'pulled'))), 
                      (aes(x=MDS1, y=MDS2, shape=soil_trt))) +
  geom_point(size=3) +
  stat_ellipse() +
  annotate('text', x=-0.8, y=0.9, label='(b) LUNA', size=8, hjust='left') +
  theme(legend.position='none')

GEMOinvaded <- ggplot(data=subset(MDScoordinates, subset=(plant_species=='GEMO' & soil_spp=='GEMO' & soil_trt %in% c('untreated', 'herbicided', 'mowed', 'pulled'))), 
                      (aes(x=MDS1, y=MDS2, shape=soil_trt))) +
  geom_point(size=3) +
  stat_ellipse() +
  annotate('text', x=-0.8, y=0.9, label='(c) GEMO', size=8, hjust='left') +
  theme(legend.position='none')

#4 panel figure
pushViewport(viewport(layout=grid.layout(1,3)))
print(ACGLinvaded, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(LUNAinvaded, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(GEMOinvaded, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export at 2000x600





