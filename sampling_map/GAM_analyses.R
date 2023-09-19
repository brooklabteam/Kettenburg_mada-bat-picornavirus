
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

#only positives that gave me sequences to work with
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
dat <- ddply(dat, .(species, month,day_of_year,year,virus,family, bat_age_class, bat_sex, specific_roost_site), summarise, N=length(any_picorna), pos=sum(any_picorna))

#get negatives and prevalence
dat$neg= datN-dat$pos
dat$prevalence <- dat$pos/dat$N

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
dat.sum$year<-factor(dat.sum$year, levels=c("2018","2019"))

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



#cbind pos neg exploration
all<-gam(cbind(pos,neg)~s(day_of_year, bs="cc"), data=dat.sum, family = binomial)
all
summary(all)
plot(all)
#not sig

species<-gam(cbind(pos,neg)~s(day_of_year, by=species, bs="cc"), data=dat.sum, family = binomial)
species
summary(species)
plot(species)
#RM and PR sig

virus<-gam(cbind(pos,neg)~s(day_of_year, by=virus, bs="cc"), data=dat.sum, family = binomial)
virus
summary(virus)
plot(virus)
#none sig

family<-gam(cbind(pos,neg)~s(day_of_year, by=family, bs="cc"), data=dat.sum, family = binomial)
family
summary(family)
plot(family)
#not sig

sex<-gam(cbind(pos,neg)~s(day_of_year, by=bat_sex, bs="cc"), data=dat.sum, family = binomial)
sex
summary(sex)
plot(sex)
#male sig

age<-gam(cbind(pos,neg)~s(day_of_year, by=bat_age_class, bs="cc"), data=dat.sum, family = binomial)
age
summary(age)
plot(age)
#not sig, adult is 0.08 so close?



#make a dataset that has individual bats listed, so we want it sorted by sampleid

homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 
dat <- read.csv(file = paste0(homewd,"/metadata/demo_sequences_picornavirales_fecesurine_meta_2018_2019.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

#clean class
unique(dat$bat_sex)

dat$bat_sex <- dat$bat_sex
dat$bat_sex[dat$bat_sex=="female"] <- "Female"
dat$bat_sex[dat$bat_sex=="male"] <- "Male"

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
#dat.list<- dlply(dat, .(sampleid))

#summarize into prevalence by species and epiwk
dat <- ddply(dat, .(sampleid, species, day_of_year,year,virus,family, bat_age_class, bat_sex, specific_roost_site), summarise, N=length(any_picorna), pos=sum(any_picorna))

#get negatives and prevalence
dat$neg= dat$N-dat$pos
dat$prevalence <- dat$pos/dat$N

#simplify= the name of "host_genus_species" 
names(dat)[names(dat)=="host_genus_species"] <- "species"

#make the data into factors
dat$species<-factor(dat$species, levels=c("Eidolon dupreanum", "Pteropus rufus", "Rousettus madagascariensis"))
dat$virus<-factor(dat$virus, levels=c("aparavirus", "bat picornavirus", "cheravirus", "cripavirus", "felisavirus",
                                              "hepatovirus","kobuvirus","kunsagivirus","mischivirus","nepovirus","picorna-like virus",
                                              "sapelovirus","sapovirus","teschovirus","tetnovirus"))
dat$family<-factor(dat$family, levels=c("picornaviridae", "caliciviridae", "secoviridae","iflaviridae","unclassified"))
dat$bat_sex<-factor(dat$bat_sex, levels=c("Male", "Female"))
dat$bat_age_class<-factor(dat$bat_age_class, levels=c("J", "A"))
dat$year<-factor(dat$year, levels=c("2018","2019"))




#individual predictors of pos neg status exploration
all<-gam(cbind(pos,neg)~s(day_of_year, bs="cc")+s(year,bs="re"), data=dat, family = binomial)
all
summary(all)
plot(all)
#not sig, no effect of year

species<-gam(prevalence~s(day_of_year,by=species, bs="cc")+s(year,bs="re"), data=dat, family = binomial)
species
summary(species)
plot(species)
#RM and PR sig, no effect of year

virus<-gam(cbind(pos,neg)~s(day_of_year,by=virus, bs="cc")+s(year,bs="re"), data=dat, family = binomial)
virus
summary(virus)
plot(virus)
#none sig, no effect of year

family<-gam(cbind(pos,neg)~s(day_of_year, by=family, bs="cc")+s(year,bs="re"), data=dat.sum, family = binomial)
family
summary(family)
plot(family)
#not sig, no effect of year

sex<-gam(cbind(pos,neg)~s(day_of_year, by=bat_sex, bs="cc")+s(year,bs="re"), data=dat.sum, family = binomial)
sex
summary(sex)
plot(sex)
#male sig, no effect of year

age<-gam(cbind(pos,neg)~s(day_of_year, by=bat_age_class,bs="cc")+s(year,bs="re"), data=dat.sum, family = binomial)
age
summary(age)
plot(age)
#not sig, no effect of year






#Make nicer looking plots


#all positives
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(all, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(all, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se

allpos<-ggplot(dat.sum)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Picornavirales prevalence over time")

allpos

#subset by species
ED<-subset(dat.sum, species=="Eidolon dupreanum")
RM<-subset(dat.sum, species=="Rousettus madagascariensis")
PR<-subset(dat.sum, species=="Pteropus rufus")

EDgam<-gam(prevalence~s(day_of_year), data=ED, family = gaussian)
RMgam<-gam(prevalence~s(day_of_year), data=RM, family = gaussian)
PRgam<-gam(prevalence~s(day_of_year), data=PR, family = gaussian)

#eidolon
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(EDgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(EDgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se

edup<-ggplot(ED)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Eidolon dupreanum")

edup


#pteropus
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(PRgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(PRgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se

pruf<-ggplot(PR)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Pteropus rufus")

pruf


#rousettus
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(RMgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(RMgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se

rmad<-ggplot(RM)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Rousettus madagascariensis")
rmad

#all three next to each other for analysis
species<-plot_grid(edup, pruf, rmad, nrow=1)
species


#subset age and sex
Adult<-subset(dat.sum, bat_age_class=="A")
Juv<-subset(dat.sum, bat_age_class=="J")
Adultgam<-gam(prevalence~s(day_of_year), data=Adult, family = gaussian)
Juvgam<-gam(prevalence~s(day_of_year), data=Juv, family = gaussian)

Fem<-subset(dat.sum, bat_sex=="Female")
Male<-subset(dat.sum, bat_sex=="Male")
Femgam<-gam(prevalence~s(day_of_year), data=Fem, family = gaussian)
Malegam<-gam(prevalence~s(day_of_year), data=Male, family = gaussian)


#females
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(Femgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Femgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se

females<-ggplot(Fem)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Females")
females


