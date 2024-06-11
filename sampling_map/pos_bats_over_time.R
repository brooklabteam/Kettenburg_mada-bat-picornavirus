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

homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 
dat <- read.csv(file = paste0(homewd,"/metadata/demo_sampleset_fecesurine_map.csv"), header = T, stringsAsFactors = F)
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
dat <- dplyr::select(dat,roost_site,bat_age_class,bat_sex, processing_date, sampling_session,
                     bat_species, sample_id, any_picorna, any_calici,any_pos, month, year)

head(dat)
unique(dat$roost_site)

#get into date form
#dat$collection_date <- as.Date(dat$collection_date,format = "%y/%m/%d")

#make the data into factors
dat$bat_species<-factor(dat$bat_species, levels=c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis"))
dat$bat_sex<-factor(dat$bat_sex, levels=c("Female","Male"))
dat$bat_age_class<-factor(dat$bat_age_class, levels=c("J", "A"))
dat$month<-factor(dat$month, levels=c("Jan", "Feb", "March", "April", "May", "June", "July", "August", "Sept",
                                      "Oct", "Nov", "Dec"))

#change picorna to numeric
dat$any_picorna <- as.numeric(dat$any_picorna)
dat$any_calici <- as.numeric(dat$any_calici)

#and make sure it is only 1 of the same sample type per each date
dat.list <- dlply(dat, .(sample_id))

#summarize into prevalence by species and epiwk
dat.sum_picorna <- ddply(dat, .(bat_species, year, bat_age_class, bat_sex, month), summarise, N=length(any_picorna), pos=sum(any_picorna))
dat.sum_calici <- ddply(dat, .(bat_species, year, bat_age_class, bat_sex, month), summarise, N=length(any_calici), pos=sum(any_calici))

#get negatives and prevalence
dat.sum_picorna$neg= dat.sum_picorna$N-dat.sum_picorna$pos
dat.sum_picorna$prevalence <- dat.sum_picorna$pos/dat.sum_picorna$N

dat.sum_calici$neg= dat.sum_calici$N-dat.sum_calici$pos
dat.sum_calici$prevalence <- dat.sum_calici$pos/dat.sum_calici$N

#and confidence intervals on the prevalence
CIs_picorna <- mapply(FUN=prop.test, x=as.list(dat.sum_picorna$pos), n=as.list(dat.sum_picorna$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)
CIs_calici <- mapply(FUN=prop.test, x=as.list(dat.sum_calici$pos), n=as.list(dat.sum_calici$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)

#and extract the upper and lower CIs
get.CIs <- function(df){
  lci = df$conf.int[1]
  uci = df$conf.int[2]
  out.df<- cbind.data.frame(lci=lci, uci=uci)
  return(out.df)
}

CIs_picorna <- lapply(CIs_picorna, get.CIs)
CIs_calici <- lapply(CIs_calici, get.CIs)

dat.sum_picorna$lci <- c(unlist(sapply(CIs_picorna, '[',1)))
dat.sum_picorna$uci <- c(unlist(sapply(CIs_picorna, '[',2)))

dat.sum_calici$lci <- c(unlist(sapply(CIs_calici, '[',1)))
dat.sum_calici$uci <- c(unlist(sapply(CIs_calici, '[',2)))

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


