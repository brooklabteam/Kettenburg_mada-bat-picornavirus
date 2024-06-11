rm(list=ls())

#packages
library(plyr) 
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(hrbrthemes)
library(reshape)
library(paletteer)
library(ggrepel)
library(hrbrthemes)
library(reshape)
library(paletteer)
library(ggh4x)
#####################################################################

homewd = "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus" 

dat <- read.csv(file = paste0(homewd,"/metadata/summary_heatmap1.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

dat$virus<-factor(dat$virus, levels=c("Bat picornavirus", "Hepatovirus", "Kobuvirus", "Kunsagivirus", 
                                      "Mischivirus", "Sapelovirus", "Teschovirus", "Sapovirus", "Aparavirus",
                                      "Cripavirus", "Cheravirus", "Nepovirus", "Felisavirus", "Picorna-like virus", "Tetnovirus"))
dat$family<-factor(dat$family, levels=c("Picornaviridae", "Caliciviridae", "Dicistroviridae", "Secoviridae", "Unclassified"))


ED<-subset(dat, species=="E. dupreanum")
PR<-subset(dat, species=="P. rufus")
RM<-subset(dat, species=="R. madagascariensis")


picorna<-subset(dat,picorna=="1")
other<-subset(dat,picorna=="0")


#facet by viral family
p9 <- ggplot(dat, aes(x=species, y=virus)) +
  geom_tile(aes(fill = cut(num_genome,breaks=0:3, labels=1:3))) +
  scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "G", direction = -1) +
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(family+type~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Count",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(face="italic"),
        legend.position = "right")
p9





p6 <- ggplot(ED, aes(x=type, y=virus)) +
  geom_tile(aes(fill = cut(num_genome,breaks=0:3, labels=1:3))) +
  scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "G", direction = -1) +
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(family~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Count",
    title="Eidolon dupreanum")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(),
        legend.position = "none")
p6



p7 <- ggplot(PR, aes(x=type, y=virus)) +
  geom_tile(aes(fill = cut(num_genome,breaks=0:3, labels=1:3))) +
  scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "G", direction = -1) +
  facet_nested(family~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "Genome type",
    y= "",
    fill="Count",
    title="Pteropus rufus")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(),
        legend.position="none")
p7


p8 <- ggplot(RM, aes(x=type, y=virus)) +
  geom_tile(aes(fill = cut(num_genome,breaks=0:3, labels=1:3))) +
  scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "G", direction = -1) +
  facet_nested(family~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Count",
    title="Rousettus madagascariensis")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(),
        legend.position="none",
        legend.margin = margin(c(0,0,0,0)))
p8


leg<-get_legend(p8)
leg<-as.ggplot(leg)

library(patchwork)
final<-p6+p7+p8
final


###Put stuff together

##Version 1

# final<- plot_grid(p6,p7,p8,labels=c("","",""),rel_widths = c(1, 1,1), rel_heights = c(1, 1,1),ncol=3, align="hv", axis="b")
# final
# 
# 
# final2<- plot_grid(p4,final,leg,labels=c("A","B",""),rel_widths = c(1, 2,0.25), rel_heights = c(2, 1),ncol=3, align="hv", axis="b")
# final2






##Try making a heatmap of the sampling data like age and percent of each bat species positive
homewd = "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus" 
dat <- read.csv(file = paste0(homewd,"/metadata/demo_picornavirales_fecesurine_meta_2018_2019.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

#clean class
unique(dat$bat_sex)

dat$bat_sex <- dat$bat_sex
dat$bat_sex[dat$bat_sex=="female"] <- "Female"
dat$bat_sex[dat$bat_sex=="male"] <- "Male"

dat$month[dat$month=="1"]<-"Jan"
dat$month[dat$month=="2"]<-"Feb"
dat$month[dat$month=="3"]<-"March"
dat$month[dat$month=="4"]<-"April"
dat$month[dat$month=="5"]<-"May"
dat$month[dat$month=="6"]<-"June"
dat$month[dat$month=="7"]<-"July"
dat$month[dat$month=="8"]<-"August"
dat$month[dat$month=="9"]<-"Sept"
dat$month[dat$month=="10"]<-"Oct"
dat$month[dat$month=="11"]<-"Nov"
dat$month[dat$month=="12"]<-"Dec"

#clean class
unique(dat$bat_age_class)

#and rank by rough age
unique(dat$young_of_year)
dat$bat_age_class <- dat$bat_age_class
dat$bat_age_class[dat$bat_age_class=="P" | dat$bat_age_class=="L"] <- "A"
dat$bat_age_class[dat$bat_age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$bat_age_class[dat$young_of_year=="yes"] <- "J"


unique(dat$any_picorna)

#select columns of interest
dat <- dplyr::select(dat,roost_site,bat_age_class,bat_sex, collection_date,
                     species, sampleid, any_picorna, month, year, day_of_year,sampleid)

head(dat)
unique(dat$roost_site)

#get into date form
dat$collection_date <- as.Date(dat$collection_date,format = "%m/%d/%y")

#change picorna to numeric
dat$any_picorna <- as.numeric(dat$any_picorna)
names(dat)[names(dat)=="bat_species"] <- "species"

#and make sure it is only 1 of the same sample type per each date
dat.list <- dlply(dat, .(sampleid))

#summarize into prevalence by species and epiwk
dat.sum <- ddply(dat, .(species, year, bat_age_class, bat_sex), summarise, N=length(any_picorna), pos=sum(any_picorna))

#get negatives and prevalence
dat.sum$neg= dat.sum$N-dat.sum$pos
dat.sum$prevalence <- dat.sum$pos/dat.sum$N

#and confidence intervals on the prevalence
CIs <- mapply(FUN=prop.test, x=as.list(dat.sum$pos), n=as.list(dat.sum$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)

#and extract the upper and lower CIs
get.CIs <- function(df){
  lci = df$conf.int[1]
  uci = df$conf.int[2]
  out.df<- cbind.data.frame(lci=lci, uci=uci)
  return(out.df)
}

CIs <- lapply(CIs, get.CIs)

dat.sum$lci <- c(unlist(sapply(CIs, '[',1)))
dat.sum$uci <- c(unlist(sapply(CIs, '[',2)))

#simplify= the name of "host_genus_species" 
names(dat.sum)[names(dat.sum)=="host_genus_species"] <- "species"

#make the data into factors
dat.sum$species<-factor(dat.sum$species, levels=c("Eidolon dupreanum", "Pteropus rufus", "Rousettus madagascariensis"))
dat.sum$bat_sex<-factor(dat.sum$bat_sex, levels=c("Male", "Female"))
dat.sum$bat_age_class<-factor(dat.sum$bat_age_class, levels=c("J", "A"))
dat.sum$month<-factor(dat.sum$month, levels=c("Jan", "Feb", "March", "April", "May", "June", "July", "August", "Sept",
                                              "Oct", "Nov", "Dec"))
#dat.sum$any_picorna<-factor(dat.sum$any_picorna, levels=c("0", "1"))

#subset some data
ED<-subset(dat.sum, species=="Eidolon dupreanum")
PR<-subset(dat.sum, species=="Pteropus rufus")
RM<-subset(dat.sum, species=="Rousettus madagascariensis")
picorna<-subset(dat.sum,any_picorna=="1")
other<-subset(dat.sum,any_picorna=="0")


#facet by year 
p10 <- ggplot(dat.sum, aes(x=species, y=bat_age_class, fill=pos)) +
  geom_tile() +
  #scale_fill_manual(values=c("lightblue1""royalblue1")) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(year+bat_sex~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Bats \n positive",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 10, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(face="italic"),
        legend.position = "right")
p10


#facet by month 
p11 <- ggplot(dat.sum, aes(x=species, y=bat_age_class, fill=pos)) +
  geom_tile() +
  #scale_fill_manual(values=c("lightblue1""royalblue1")) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(month+bat_sex~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Bats \n positive",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(face="italic"),
        legend.position = "right")
p11
