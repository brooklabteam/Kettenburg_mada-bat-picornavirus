rm(list=ls())

#packages
library(sf)
library(mapplots)
library(scatterpie)
library(maptools)
library(plyr) 
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ggspatial)
library(ggrepel)
library(ggridges)

#####################################################################

homewd = "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus" 

dat <- read.csv(file = paste0(homewd,"/metadata/demo_picornavirales_fecesurine_meta_ridgeline_2018_2019.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

unique(dat$any_picorna)


#and rank by rough age, by total positives
unique(dat$young_of_year)
dat$bat_age_class <- dat$bat_age_class
dat$bat_age_class[dat$bat_age_class=="P" | dat$bat_age_class=="L"] <- "A"
dat$bat_age_class[dat$bat_age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$bat_age_class[dat$young_of_year=="yes"] <- "J"

#get into date form
dat$collection_date <- as.Date(dat$collection_date,format = "%m/%d/%y")


#check sites
unique(dat$roost_site)


#change picorna to numeric
dat$any_picorna <- as.numeric(dat$any_picorna)

names(dat)[names(dat)=="bat_species"] <- "species"

#and make sure it is only 1 of the same sample type per each date
dat.list <- dlply(dat, .(sampleid))

#summarize into prevalence by species and epiwk
dat.sum <- ddply(dat, .(species,month, year), summarise, N=length(any_picorna), pos=sum(any_picorna))

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

dat.sum$year<-factor(dat.sum$year)
dat.sum$year<-factor(dat.sum$year, levels = c("2018","2019"))

#here's a vector assigning colors to each species
colz = c("Eidolon dupreanum"="coral2", "Pteropus rufus" = "cornflowerblue", "Rousettus madagascariensis" = "darkgoldenrod1" )
colz2 = c("2018"="cornflowerblue", "2019"="darkgoldenrod1")


##Make two plots each of 2018 and 2019 just so we can see the prevalence on the y scale


##Version using lines and points
p1 <- ggplot(dat.sum, aes(x=month,y=prevalence, fill=year))+
  #geom_area(aes(x=day_of_year, y= prevalence, color=year))+
  geom_point(aes(x=month, y= prevalence, color=year, size=N)) +
  geom_errorbar(aes(x=month, ymin=lci, ymax=uci, color=year), size=.1) +
  geom_smooth(aes(x=month, color=year,fill=year), se=FALSE)+
  #geom_line(aes(x=month, y= prevalence, color=year))+
  # annotate("rect", xmin = 152, xmax = 244, ymin = 0, ymax = 1,
  #          alpha = .3)+
  # geom_vline(xintercept=c(152,244), linetype="dotted")+
  scale_x_continuous(limits = c(1,12), expand=c(0,0))+
  #facet_wrap(~year, dir = "v")+
  theme_linedraw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=colz2)+
  scale_color_manual(values=colz2) +
  labs(x="Month", y="Prevalence", title = expression(~italic(Picornavirales) ~"positive bats over time"))
p1


