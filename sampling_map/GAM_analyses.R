
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
library(hrbrthemes)
library(reshape)
library(paletteer)
library(ggh4x)
library(cowplot)
library(ggplotify)
library(mgcv)


homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus"

#all positives for any picorna contig hit
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 
dat <- read.csv(file = paste0(homewd,"/metadata/demo_sequences_picornavirales_fecesurine_meta_2018_2019.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

#clean class
unique(dat$bat_sex)

dat$bat_sex <- dat$bat_sex
dat$bat_sex[dat$bat_sex=="female"] <- "Female"
dat$bat_sex[dat$bat_sex=="male"] <- "Male"

# dat$month[dat$month=="1"]<-"Jan"
# dat$month[dat$month=="2"]<-"Feb"
# dat$month[dat$month=="3"]<-"March"
# dat$month[dat$month=="4"]<-"April"
# dat$month[dat$month=="5"]<-"May"
# dat$month[dat$month=="6"]<-"June"
# dat$month[dat$month=="7"]<-"July"
# dat$month[dat$month=="8"]<-"August"
# dat$month[dat$month=="9"]<-"Sept"
# dat$month[dat$month=="10"]<-"Oct"
# dat$month[dat$month=="11"]<-"Nov"
# dat$month[dat$month=="12"]<-"Dec"

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
                     species, sampleid, any_picorna,virus, type, family,
                     month, year, day_of_year, specific_roost_site)

head(dat)
unique(dat$roost_site)

#get into date form
dat$collection_date <- as.Date(dat$collection_date,format = "%m/%d/%y")

#change picorna to numeric
dat$any_picorna <- as.numeric(dat$any_picorna)
names(dat)[names(dat)=="bat_species"] <- "species"

#and make sure it is only 1 of the same sample type per each date
#dat.list <- dlply(dat, .(sampleid))

#summarize into prevalence by species and epiwk
dat.sum <- ddply(dat, .(species, month,day_of_year,year,virus,family, bat_age_class, bat_sex, specific_roost_site), summarise, N=length(any_picorna), pos=sum(any_picorna))

#get negatives and prevalence
dat.sum$neg= dat.sum$N-dat.sum$pos
dat.sum$prevalence <- dat.sum$pos/dat.sum$N

#simplify= the name of "host_genus_species" 
names(dat.sum)[names(dat.sum)=="host_genus_species"] <- "species"

#make the data into factors
dat.sum$species<-factor(dat.sum$species, levels=c("Eidolon dupreanum", "Pteropus rufus", "Rousettus madagascariensis"))
dat.sum$virus<-factor(dat.sum$virus, levels=c("aparavirus", "bat picornavirus", "cheravirus", "cripavirus", "felisavirus",
                                              "hepatovirus","kobuvirus","kunsagivirus","mischivirus","nepovirus","picorna-like virus",
                                              "sapelovirus","sapovirus","teschovirus","tetnovirus"))
dat.sum$family<-factor(dat.sum$family, levels=c("picornaviridae", "caliciviridae", "secoviridae","iflaviridae","unclassified"))
dat.sum$bat_sex<-factor(dat.sum$bat_sex, levels=c("Male", "Female"))
dat.sum$bat_age_class<-factor(dat.sum$bat_age_class, levels=c("J", "A"))


#dat.sum$any_picorna<-factor(dat.sum$any_picorna, levels=c("0", "1"))


#prevalence by day of year for all species
plot(prevalence~day_of_year, data=dat.sum, type="l")
all<-gam(prevalence~s(day_of_year), data=dat.sum, family = gaussian)
all
summary(all)

plot(all, shade=TRUE)
#not significant

#by bat species
bat<-gam(prevalence~s(day_of_year, by=species), data=dat.sum, family = gaussian)
bat
summary(bat)

plot(bat, shade=TRUE)
#RM is significant

#by age and sex
age<-gam(prevalence~s(day_of_year, by=bat_age_class), data=dat.sum, family = gaussian)
age
summary(age)

plot(age, shade=TRUE)

sex<-gam(prevalence~s(day_of_year, by=bat_sex), data=dat.sum, family = gaussian)
sex
summary(sex)

plot(sex, shade=TRUE)

#age not significant, sex not significant

#by virus genus and family
virus<-gam(prevalence~s(day_of_year, by=virus), data=dat.sum, family = gaussian)
virus
summary(virus)

plot(virus, shade=TRUE)

family<-gam(prevalence~s(day_of_year, by=family), data=dat.sum, family = gaussian)
family
summary(family)

plot(family, shade=TRUE)
#picorna-like virus and teschovirus significant, unclassified significant

