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
dat <- read.csv(file = paste0(homewd, "fullgenome_coverage.csv"), header = T, stringsAsFactors = F) #load data
dat$virus<-factor(dat$virus, levels=c("Icavirus", "Kobuvirus", "Kunsagivirus", "Nepovirus", "Picornavirus", "Sapelovirus",
                                      "Sapovirus", "Teschovirus")) ##pick order to display
dat$virus<-factor(dat$species, levels=c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis")) ##pick order to display

##subset the data into different viruses
icavirus<-subset(dat, virus=="Icavirus")
kobuvirus<-subset(dat, virus=="Kobuvirus")
kunsagivirus<-subset(dat, virus=="Kunsagivirus")
nepovirus<-subset(dat, virus=="Nepovirus")
nepovirusRNA1<-subset(nepovirus, accession=="OQ818326: Nepovirus RNA 1")
nepovirusRNA2<-subset(nepovirus, accession=="OQ818327: Nepovirus RNA 2")
picornavirus<-subset(dat, virus=="Picornavirus")
picornavirus1<-subset(picornavirus, accession=="OQ818325: Picornavirus")
picornavirus2<-subset(picornavirus, accession=="OQ818328: Picornavirus")
sapelovirus<-subset(dat, virus=="Sapelovirus")
sapelovirus1<-subset(sapelovirus, accession=="OQ818329: Sapelovirus")
sapelovirus2<-subset(sapelovirus, accession=="OQ818320: Sapelovirus")
sapelovirus3<-subset(sapelovirus, accession=="OQ818321: Sapelovirus")
sapovirus<-subset(dat, virus=="Sapovirus")
teschovirus<-subset(dat, virus=="Teschovirus")
teschovirus1<-subset(teschovirus, accession=="OQ818323: Teschovirus")
teschovirus2<-subset(teschovirus, accession=="OQ818324: Teschovirus")


dat <- read.csv(file = paste0(homewd, "allgenome_coverage.csv"), header = T, stringsAsFactors = F) #load data
dat$virus<-factor(dat$virus, levels=c("Icavirus", "Kobuvirus", "Kunsagivirus", "Nepovirus", "Picornavirus", "Sapelovirus",
                                      "Sapovirus", "Teschovirus")) ##pick order to display
dat$species<-factor(dat$species, levels=c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis")) ##pick order to display
dat$type<-factor(dat$type, levels=c("full", "partial"))

#subset into diff bat species
pr<-subset(dat, species=="Pteropus rufus")
ed<-subset(dat, species=="Eidolon dupreanum")
rm<-subset(dat, species=="Rousettus madagascariensis")


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


colz=c("partial"="deepskyblue3", "full"="coral2")
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


colz=c("partial"="deepskyblue3", "full"="coral2")
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



all<-plot_grid(pr_all,ed_all, rm_all, ncol=3, labels = "AUTO", align="hv", axis="l")
all






##If you want to plot individual viruses and stitch together that way


#plot data by rPM
title<-expression(paste(italic("Pteropus rufus mischivirus "), "OQ818316"))
mischi <- ggplot(icavirus, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
mischi


title<-expression(paste(italic("Eidolon dupreanum kobuvirus "), "OQ818322"))
kobu <- ggplot(kobuvirus, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
kobu


title<-expression(paste(italic("Eidolon dupreanum kunsagivirus "), "OQ818317"))
kun <- ggplot(kunsagivirus, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
kun


title<-expression(paste(italic("Rousettus madagascariensis nepovirus")))
nepo <- ggplot(nepovirus, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue4", "deepskyblue"))+
  theme_linedraw() +
  theme(legend.position = c(0.8,0.7))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  theme(legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'))+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
nepo


title<-expression(paste(italic("Rousettus madagascariensis picornavirus "), "OQ818325"))
picorna1 <- ggplot(picornavirus1, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
picorna1


title<-expression(paste(italic("Rousettus madagascariensis picornavirus "), "OQ818328"))
picorna2 <- ggplot(picornavirus2, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
picorna2


title<-expression(paste(italic("Eidolon dupreanum sapelovirus "), "OQ818320"))
sapelo1 <- ggplot(sapelovirus1, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapelo1


title<-expression(paste(italic("Eidolon dupreanum sapelovirus "), "OQ818321"))
sapelo2 <- ggplot(sapelovirus2, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapelo2


title<-expression(paste(italic("Rousettus madagascariensis sapelovirus "), "OQ818329"))
sapelo3 <- ggplot(sapelovirus3, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapelo3


title<-expression(paste(italic("Eidolon dupreanum sapovirus "), "OQ818319"))
sapo <- ggplot(sapovirus, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapo


title<-expression(paste(italic("Eidolon dupreanum teschovirus "), "OQ818318"))
tescho1 <- ggplot(teschovirus1, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
tescho1


title<-expression(paste(italic("Rousettus madagascariensis teschovirus "), "OQ818323"))
tescho2 <- ggplot(teschovirus2, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(trans="log10",expand = c(0,0), limits = c(1,10000))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million) log10 scale")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
tescho2


title<-expression(paste(italic("Rousettus madagascariensis teschovirus "), "OQ818324"))
tescho3 <- ggplot(teschovirus1, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("deepskyblue3"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
tescho3






##next do partial genomes
dat <- read.csv(file = paste0(homewd, "partialgenome_coverage.csv"), header = T, stringsAsFactors = F) #load data
dat$virus<-factor(dat$virus, levels=c("Bee virus", "Cheravirus", "Cripavirus", "Felisavirus", "Hepatovirus", "Picorna-like virus",
                                      "Picornavirus", "Sapelovirus","Sapovirus", "Tetnovirus")) #pick order to display


##subset the data into different viruses
bee<-subset(dat, virus=="Bee virus")
chera<-subset(dat, virus=="Cheravirus")
cripa<-subset(dat, virus=="Cripavirus")
felisa<-subset(dat, virus=="Felisavirus")
felisa1<-subset(felisa, accession=="OQ818335: Felisavirus")
felisa2<-subset(felisa, accession=="OQ818341: Felisavirus")
hep<-subset(dat, virus=="Hepatovirus")
picornalike<-subset(dat, virus=="Picorna-like virus")
picornalike1<-subset(picornalike, accession=="OQ818333: Picorna-like virus")
picornalike2<-subset(picornalike, accession=="OQ818334: Picorna-like virus")
picornalike3<-subset(picornalike, accession=="OQ818336: Picorna-like virus")
picornalike4<-subset(picornalike, accession=="OQ818339: Picorna-like virus")
picorna_partial<-subset(dat, virus=="Picornavirus")
sapelo_partial<-subset(dat, virus=="Sapelovirus")
sapelovirus1partial<-subset(sapelo_partial, accession=="OQ818342: Sapelovirus: Sapelovirus")
sapelovirus2partial<-subset(sapelo_partial, accession=="OQ818343: Sapelovirus")
sapelovirus3partial<-subset(sapelo_partial, accession=="OQ818344: Sapelovirus")
sapovirus_partial<-subset(dat, virus=="Sapovirus")
sapovirus1partial<-subset(sapovirus_partial, accession=="OQ818340: Sapovirus")
sapovirus2partial<-subset(sapovirus_partial, accession=="OQ818345: Sapovirus")
sapovirus3partial<-subset(sapovirus_partial, accession=="OQ818347: Sapovirus")
sapovirus4partial<-subset(sapovirus_partial, accession=="OQ818348: Sapovirus")
tetnovirus<-subset(dat, virus=="Tetnovirus")



#plot data by rPM
title<-expression(paste(italic("Pteropus rufus bee virus "), "OQ818332"))
beevirus <- ggplot(bee, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
beevirus


title<-expression(paste(italic("Pteropus rufus cheravirus "), "OQ818330"))
cheravirus <- ggplot(chera, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
cheravirus


title<-expression(paste(italic("Eidolon dupreanum cripavirus "), "OQ818338"))
cripavirus <- ggplot(cripa, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
cripavirus


title<-expression(paste(italic("Pteropus rufus felisavirus "), "OQ818335"))
felisavirus1 <- ggplot(felisa1, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
felisavirus1


title<-expression(paste(italic("Eidolon dupreanum felisavirus "), "OQ818341"))
felisavirus2 <- ggplot(felisa2, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
felisavirus2


title<-expression(paste(italic("Eidolon dupreanum hepatovirus "), "OQ818337"))
hepatovirus <- ggplot(hep, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
hepatovirus


title<-expression(paste(italic("Pteropus rufus picorna-like virus "), "OQ818333"))
picornalike1 <- ggplot(picornalike1, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
picornalike1


title<-expression(paste(italic("Pteropus rufus picorna-like virus "), "OQ818334"))
picornalike2 <- ggplot(picornalike2, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
picornalike2

title<-expression(paste(italic("Eidolon dupreanum picorna-like virus "), "OQ818336"))
picornalike3 <- ggplot(picornalike3, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
picornalike3


title<-expression(paste(italic("Eidolon dupreanum picorna-like virus "), "OQ818339"))
picornalike4 <- ggplot(picornalike4, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
picornalike4



title<-expression(paste(italic("Rousettus madagascariensis picornavirus "), "OQ818346"))
picornapartial <- ggplot(picorna_partial, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
picornapartial


title<-expression(paste(italic("Eidolon dupreanum sapelovirus "), "OQ818342"))
sapelopartial1 <- ggplot(sapelovirus1partial, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapelopartial1


title<-expression(paste(italic("Eidolon dupreanum sapelovirus "), "OQ818343"))
sapelopartial2 <- ggplot(sapelovirus2partial, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapelopartial2


title<-expression(paste(italic("Eidolon dupreanum sapelovirus "), "OQ818344"))
sapelopartial3 <- ggplot(sapelovirus3partial, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapelopartial3


title<-expression(paste(italic("Eidolon dupreanum sapovirus "), "OQ818340"))
sapopartial1 <- ggplot(sapovirus1partial, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapopartial1


title<-expression(paste(italic("Rousettus madagascariensis sapovirus "), "OQ818345"))
sapopartial2 <- ggplot(sapovirus2partial, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapopartial2


title<-expression(paste(italic("Rousettus madagascariensis sapovirus "), "OQ818347"))
sapopartial3 <- ggplot(sapovirus3partial, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapopartial3


title<-expression(paste(italic("Rousettus madagascariensis sapovirus "), "OQ818348"))
sapopartial4 <- ggplot(sapovirus4partial, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sapopartial4


title<-expression(paste(italic("Pteropus rufus tetnovirus "), "OQ818331"))
tetno <- ggplot(tetnovirus, aes(x=Position, y=rPM, fill=accession))+ 
  geom_area(linewidth=0.5) +
  scale_fill_manual(values=c("coral2"))+
  theme_linedraw() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle(title)+
  xlab("Genome position") + ylab("rPM (reads per million)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
tetno




#sort by bat species
library(patchwork)
library(cowplot)
library(ggplotify)


p_ruf<-mischi+cheravirus+tetno+beevirus+picornalike1+picornalike2+picornalike3+felisavirus1+plot_layout(nrow=1) #8 viruses
p_ruf

e_dup<-kun+tescho1+sapo+sapelo1+sapelo2+kobu+hepatovirus+cripavirus+picornalike4+sapopartial1+
  felisavirus2+sapelopartial1+sapelopartial2+sapelopartial3+plot_layout(nrow=4, ncol=4) #14 viruses
e_dup

r_mad<-tescho2+tescho3+picorna1+picorna2+nepo+sapelo3+picornapartial+sapopartial2+sapopartial3+sapopartial4+
  plot_layout(nrow=3, ncol=4) #10 viruses
r_mad

all<-plot_grid(e_dup,r_mad, ncol=8, nrow=5, rel_widths = c(1,1), rel_heights = c(1,1), align = "hv", axis="b")
all

final<-plot_grid(p_ruf,all, nrow=2, rel_widths = c(5,1), rel_heights = c(1,5), align = "hv", axis="b")
final

p_ruf_full<-mischi+plot_spacer()+plot_spacer()+plot_spacer()+plot_spacer()+plot_spacer()+plot_layout(nrow=1)
p_ruf_full
p_ruf_full<-as.ggplot(p_ruf_full)

r_mada_full<-tescho2+tescho3+picorna1+picorna2+nepo+sapelo3+plot_layout(nrow=1)
r_mada_full
r_mada_full<-as.ggplot(r_mada_full)

e_dup_full<-kun+tescho1+sapo+sapelo1+sapelo2+kobu+plot_layout(nrow=1)
e_dup_full
e_dup_full<-as.ggplot(e_dup_full)


full<-e_dup_full/r_mada_full/p_ruf_full+plot_layout(nrow=3, ncol = 1, heights=c(1,1,1), widths=c(1,1,0.5))+plot_annotation(tag_levels = "A")
full
full<-as.ggplot(full)


seq1<-cowplot::plot_grid(e_dup_full,r_mada_full,p_ruf_full,nrow = 3, ncol=1, rel_heights =c(1,1,1), rel_widths = c(1,1,1), labels = "AUTO", align="hv", axis="l", scale = c(1,1,1))
seq1
