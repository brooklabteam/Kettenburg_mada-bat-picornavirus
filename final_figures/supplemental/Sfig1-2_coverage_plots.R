#This is to plot a visual supplemental fig of the coverage plots for the genomes

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

homewd="/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/read_mapping/"
setwd(paste0(homewd))

dat <- read.csv(file= "all_genome_coverage.csv", header = T, stringsAsFactors = F) #load data
dat$virus<-factor(dat$virus, levels=c("Mischivirus", "Kobuvirus", "Kunsagivirus", "Nepovirus", "Roupivirus", "Sapelovirus",
                                      "Sapovirus", "Teschovirus")) ##pick order to display
dat$species<-factor(dat$species, levels=c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis")) ##pick order to display
dat$type<-factor(dat$type, levels=c("full", "partial"))

#subset into diff bat species and bu full and partial genomes
pr<-subset(dat, species=="Pteropus rufus")
ed<-subset(dat, species=="Eidolon dupreanum")
rm<-subset(dat, species=="Rousettus madagascariensis")

pr_full<-subset(pr, type=="full")
ed_full<-subset(ed, type=="full")
rm_full<-subset(rm, type=="full")

ed_partial<-subset(ed, type=="partial")
rm_partial<-subset(rm, type=="partial")

colz=c("partial"="deepskyblue3", "full"="coral2")

#just the full genomes

#sub out reads per million and coverage, then make both figures
title<-expression(paste(italic("Pteropus rufus")))
#pr_f<-ggplot(pr_full, aes(x=Position, y=rPM, fill=type))+ 
pr_f<-ggplot(pr_full, aes(x=Position, y=Coverage, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  #xlab("Genome position") + ylab("rPM (reads per million)")+
  xlab("Genome position") + ylab("Coverage (reads)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
pr_f


title<-expression(paste(italic("Eidolon dupreanum")))
#ed_f<-ggplot(ed_full, aes(x=Position, y=rPM, fill=type))+ 
ed_f<-ggplot(ed_full, aes(x=Position, y=Coverage, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  #xlab("Genome position") + ylab("rPM (reads per million)")+
  xlab("Genome position") + ylab("Coverage (reads)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
ed_f


title<-expression(paste(italic("Rousettus madagascariensis")))
#rm_f<-ggplot(rm_full, aes(x=Position, y=rPM, fill=type))+ 
rm_f<-ggplot(rm_full, aes(x=Position, y=Coverage, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("deepskyblue"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  #xlab("Genome position") + ylab("rPM (reads per million)")+
  xlab("Genome position") + ylab("Coverage (reads)")+
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


#just the partial genomes
title<-expression(paste(italic("Eidolon dupreanum")))
#ed_p<-ggplot(ed_partial, aes(x=Position, y=rPM, fill=type))+ 
ed_p<-ggplot(ed_partial, aes(x=Position, y=Coverage, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  #xlab("Genome position") + ylab("rPM (reads per million)")+
  xlab("Genome position") + ylab("Coverage (reads)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
ed_p


title<-expression(paste(italic("Rousettus madagascariensis")))
#rm_p<-ggplot(rm_partial, aes(x=Position, y=rPM, fill=type))+ 
rm_p<-ggplot(rm_partial, aes(x=Position, y=Coverage, fill=type))+ 
  geom_area(linewidth=0.5) +
  facet_nested(accession~.,
               scales="free", nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  ggtitle(title)+
  scale_y_continuous(position="left")+
  #xlab("Genome position") + ylab("rPM (reads per million)")+
  xlab("Genome position") + ylab("Coverage (reads)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(0,0,0,1)))
rm_p


partial<-plot_grid(ed_p,rm_p, ncol=2, labels =c("A","B"),label_size=23, align="hv", axis="l")
partial


# save figs
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 


#ggsave(file = paste0(homewd, "/figures/full_coverage_rpm.pdf"),
ggsave(file = paste0(homewd, "/final_figures/supplemental/SfigX_full_coverage_plot.pdf"),
       plot= full,
       units="mm",  
       width=200, 
       height=110, 
       scale=2, 
       dpi=500)

ggsave(file = paste0(homewd, "/final_figures/supplemental/SfigX_partial_coverage_plot.pdf"),
       plot= partial,
       units="mm",  
       width=220, 
       height=220, 
       scale=2, 
       dpi=500)

