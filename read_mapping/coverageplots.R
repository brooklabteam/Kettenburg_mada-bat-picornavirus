rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gggenes)
library(devtools)
library(gridExtra)
library(cowplot)
library(patchwork)
library(ggh4x)

## now let's plot coverage 

homewd="/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/read_mapping/"
setwd(paste0(homewd))

##first do full genomes
dat <- read.csv(file = "fullgenome_coverage.csv", header = T, stringsAsFactors = F) #load data
dat$virus<-factor(dat$virus, levels=c("Mischivirus", "Kobuvirus", "Kunsagivirus", "Nepovirus", "Roupivirus", "Sapelovirus",
                                      "Saplivirus", "Teschovirus")) ##pick order to display
dat$virus<-factor(dat$species, levels=c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis")) ##pick order to display

##subset the data into different viruses
Mischivirus<-subset(dat, virus=="Mischivirus")
kobuvirus<-subset(dat, virus=="Kobuvirus")
kunsagivirus<-subset(dat, virus=="Kunsagivirus")
nepovirus<-subset(dat, virus=="Nepovirus")
nepovirusRNA1<-subset(nepovirus, accession=="OQ818326: Nepovirus RNA 1")
nepovirusRNA2<-subset(nepovirus, accession=="OQ818327: Nepovirus RNA 2")
Roupivirus<-subset(dat, virus=="Roupivirus")
Roupivirus1<-subset(Roupivirus, accession=="OQ818325: Roupivirus")
Roupivirus2<-subset(Roupivirus, accession=="OQ818328: Roupivirus")
sapelovirus<-subset(dat, virus=="Sapelovirus")
sapelovirus1<-subset(sapelovirus, accession=="OQ818329: Sapelovirus")
sapelovirus2<-subset(sapelovirus, accession=="OQ818320: Sapelovirus")
sapelovirus3<-subset(sapelovirus, accession=="OQ818321: Sapelovirus")
Saplivirus<-subset(dat, virus=="Saplivirus")
teschovirus<-subset(dat, virus=="Teschovirus")
teschovirus1<-subset(teschovirus, accession=="OQ818323: Teschovirus")
teschovirus2<-subset(teschovirus, accession=="OQ818324: Teschovirus")


dat <- read.csv(file= "allgenome_coverage.csv", header = T, stringsAsFactors = F) #load data
dat$virus<-factor(dat$virus, levels=c("Mischivirus", "Kobuvirus", "Kunsagivirus", "Nepovirus", "Roupivirus", "Sapelovirus",
                                      "Saplivirus", "Teschovirus")) ##pick order to display
dat$species<-factor(dat$species, levels=c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis")) ##pick order to display
dat$type<-factor(dat$type, levels=c("full", "partial"))

#subset into diff bat species and bu full and partial genomes
pr<-subset(dat, species=="Pteropus rufus")
ed<-subset(dat, species=="Eidolon dupreanum")
rm<-subset(dat, species=="Rousettus madagascariensis")

pr_full<-subset(pr, type=="full")
ed_full<-subset(ed, type=="full")
rm_full<-subset(rm, type=="full")

pr_partial<-subset(pr, type=="partial")
ed_partial<-subset(ed, type=="partial")
rm_partial<-subset(rm, type=="partial")

colz=c("partial"="deepskyblue3", "full"="coral2")

title<-expression(paste(italic("Pteropus rufus")))
pr_all<-ggplot(pr, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue","coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
pr_all


title<-expression(paste(italic("Eidolon dupreanum")))
ed_all<-ggplot(ed, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue","coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
ed_all


title<-expression(paste(italic("Rousettus madagascariensis")))
rm_all<-ggplot(rm, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue","coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
rm_all

#everything together
all<-plot_grid(pr_all,ed_all, rm_all, ncol=3, labels = "AUTO", align="hv", axis="l")
all

#export 15x20 inch portrait PDF


#just the full genomes
title<-expression(paste(italic("Pteropus rufus")))
pr_f<-ggplot(pr_full, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
pr_f


title<-expression(paste(italic("Eidolon dupreanum")))
ed_f<-ggplot(ed_full, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
ed_f


title<-expression(paste(italic("Rousettus madagascariensis")))
rm_f<-ggplot(rm_full, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
rm_f

library(patchwork)

ed_pf<-plot_grid(pr_f,ed_f, ncol=1,labels="AUTO", align="hv", axis="l")
ed_pf

ed_pf<-pr_f/ed_f+plot_layout(ncol=1,  heights = c(0.15, 1))+plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 23))&theme(plot.tag = element_text(face = 'bold'))
ed_pf

ed_pf<-as.ggplot(ed_pf)


full<-plot_grid(ed_pf, rm_f, ncol=2, labels = c("","C"),label_size=23, align="hv", axis="l")
full

#export pdf 18x13 inch PDF landscape

#just the partial genomes
title<-expression(paste(italic("Pteropus rufus")))
pr_p<-ggplot(pr_partial, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
pr_p


title<-expression(paste(italic("Eidolon dupreanum")))
ed_p<-ggplot(ed_partial, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
ed_p


title<-expression(paste(italic("Rousettus madagascariensis")))
rm_p<-ggplot(rm_partial, aes(x=Position, y=rPM, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  xlab("Genome position") + ylab("rPM (reads per million)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
rm_p


partial<-plot_grid(pr_p, ed_p,rm_p, ncol=3, labels =c("A","B","C"),label_size=23, align="hv", axis="l")
partial

#export 15x13inch PDF landscape

