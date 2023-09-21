
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


#only positives that gave me sequences to work with
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 
dat <- read.csv(file = paste0(homewd,"/metadata/demo_sequences_picornavirales_fecesurine_meta_2018_2019_byday.csv"), header = T, stringsAsFactors = F)
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
#dat.list <- dlply(dat, .(sampleid))

#summarize into prevalence by species and epiwk
dat <- ddply(dat, .(species,day_of_year,year,virus,family, bat_age_class, bat_sex, specific_roost_site), summarise, N=length(any_picorna), pos=sum(any_picorna))

#get negatives and prevalence
dat$neg= dat$N-dat$pos
dat$prevalence <- dat$pos/dat$N

#simplify= the name of "host_genus_species" 
names(dat)[names(dat)=="host_genus_species"] <- "species"
dat$virus  #some are blank, so convert those to NA
dat$virus[dat$virus==""] <- NA
dat$family  #some are blank, so convert those to NA
dat$family[dat$virus==""] <- NA

#make the data into factors
dat$species<-factor(dat$species, levels=c("Eidolon dupreanum", "Pteropus rufus", "Rousettus madagascariensis"))
dat$virus<-factor(dat$virus, levels=c("aparavirus", "bat picornavirus", "cheravirus", "cripavirus", "felisavirus",
                                              "hepatovirus","kobuvirus","kunsagivirus","mischivirus","nepovirus","picorna-like virus",
                                              "sapelovirus","sapovirus","teschovirus","tetnovirus"))
dat$family<-factor(dat$family, levels=c("picornaviridae", "caliciviridae", "secoviridae","iflaviridae","unclassified"))
dat$bat_sex<-factor(dat$bat_sex, levels=c("Male", "Female"))
dat$bat_age_class<-factor(dat$bat_age_class, levels=c("J", "A"))
dat$year<-factor(dat$year, levels=c("2018","2019"))

#dat$any_picorna<-factor(dat$any_picorna, levels=c("0", "1"))


#prevalence by day of year for all species
plot(prevalence~day_of_year, data=dat, type="l")
all<-gam(prevalence~s(day_of_year), data=dat, family = gaussian)
all
summary(all)
plot(all, shade=TRUE)
#not significant

#by bat species
bat<-gam(prevalence~s(day_of_year, by=species), data=dat, family = gaussian)
bat
summary(bat)
plot(bat, shade=TRUE)
#RM is significant

#by age and sex
age<-gam(prevalence~s(day_of_year, by=bat_age_class), data=dat, family = gaussian)
age
summary(age)
plot(age, shade=TRUE)

sex<-gam(prevalence~s(day_of_year, by=bat_sex), data=dat, family = gaussian)
sex
summary(sex)
plot(sex, shade=TRUE)

#age not significant, sex not significant

#by virus genus and family
virus<-gam(prevalence~s(day_of_year, by=virus), data=dat, family = gaussian)
virus
summary(virus)

plot(virus, shade=TRUE)

family<-gam(prevalence~s(day_of_year, by=family), data=dat, family = gaussian)
family
summary(family)
plot(family, shade=TRUE)
#teschovirus and picorna-like virus significant, unclassified significant



#cbind pos neg exploration
all<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=dat, family = binomial)
all
summary(all)
plot(all)
#not sig

species<-gam(cbind(pos,neg)~s(day_of_year, by=species, bs="cc", k=7), data=dat, family = binomial)
species
summary(species)
plot(species)
#RM and PR sig

virus<-gam(cbind(pos,neg)~s(day_of_year, by=virus, bs="cc", k=7), data=dat, family = binomial)
virus
summary(virus)
plot(virus)
#none sig

family<-gam(cbind(pos,neg)~s(day_of_year, by=family, bs="cc", k=7), data=dat, family = binomial)
family
summary(family)
plot(family)
#not sig

sex<-gam(cbind(pos,neg)~s(day_of_year, by=bat_sex, bs="cc", k=7), data=dat, family = binomial)
sex
summary(sex)
plot(sex)
#male sig

age<-gam(cbind(pos,neg)~s(day_of_year, by=bat_age_class, bs="cc", k=7), data=dat, family = binomial)
age
summary(age)
plot(age)
#not sig



#Make nicer looking plots for "by day"


#all positives
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(all, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(all, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

allpos<-ggplot(dat)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Picornavirales prevalence over time")

allpos

#subset by species
ED<-subset(dat, species=="Eidolon dupreanum")
RM<-subset(dat, species=="Rousettus madagascariensis")
PR<-subset(dat, species=="Pteropus rufus")

EDgam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=ED, family = binomial)
RMgam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=RM, family = binomial)
PRgam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=PR, family = binomial)

#eidolon
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(EDgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(EDgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

edup<-ggplot(ED)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
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
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

pruf<-ggplot(PR)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
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
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

rmad<-ggplot(RM)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
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
Adult<-subset(dat, bat_age_class=="A")
Juv<-subset(dat, bat_age_class=="J")
Adultgam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=Adult, family = binomial)
Juvgam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=Juv, family = binomial)

Fem<-subset(dat, bat_sex=="Female")
Male<-subset(dat, bat_sex=="Male")
Femgam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=Fem, family = binomial)
Malegam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=Male, family = binomial)


#females
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(Femgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Femgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

females<-ggplot(Fem)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Females")
females


#males
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(Malegam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Malegam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

males<-ggplot(Male)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Males")
males

#adults
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(Adultgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Adultgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

adults<-ggplot(Adult)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Adults")
adults


#juveniles
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(Juvgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Juvgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

juveniles<-ggplot(Juv)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Juveniles")
juveniles

#plot next to each other for analysis
sexage<-plot_grid(females, males,adults, juveniles, nrow=2)
sexage


#subset by viral family
picorna<-subset(dat, family=="picornaviridae")
cali<-subset(dat, family=="caliciviridae")
seco<-subset(dat, family=="secoviridae")
ifla<-subset(dat, family=="iflaviridae")
un<-subset(dat, family=="unclassified")

picornagam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7), data=picorna, family = binomial)
caligam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=5), data=cali, family = binomial)
ungam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=4), data=un, family = binomial)

#picornaviridae
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(picornagam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(picornagam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

picornaviridae<-ggplot(picorna)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Picornaviridae")
picornaviridae

#caliciviridae
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(caligam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(caligam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

caliciviridae<-ggplot(cali)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Caliciviridae")
caliciviridae

#unclassified
pred.df=cbind.data.frame(day_of_year=1:365)
pred.df$predicted_count<-predict.gam(ungam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(ungam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

unclassified<-ggplot(un)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Unclassified family")
unclassified

#plot together
families<-plot_grid(picornaviridae,caliciviridae,unclassified, nrow=1)
families






#make a dataset that has individual bats listed, so we want it sorted by sampleid
rm(list=ls())
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 
dat2 <- read.csv(file = paste0(homewd,"/metadata/demo_sequences_picornavirales_fecesurine_meta_2018_2019_byday.csv"), header = T, stringsAsFactors = F)
head(dat2)
names(dat2)

#clean class
unique(dat2$bat_sex)

dat2$bat_sex <- dat2$bat_sex
dat2$bat_sex[dat2$bat_sex=="female"] <- "Female"
dat2$bat_sex[dat2$bat_sex=="male"] <- "Male"

#clean class
unique(dat2$bat_age_class)

#and rank by rough age
unique(dat2$young_of_year)
dat2$bat_age_class <- dat2$bat_age_class
dat2$bat_age_class[dat2$bat_age_class=="P" | dat2$bat_age_class=="L"] <- "A"
dat2$bat_age_class[dat2$bat_age_class=="NL" | dat2$young_of_year=="no"] <- "A"
dat2$bat_age_class[dat2$young_of_year=="yes"] <- "J"

unique(dat2$any_picorna)

#select columns of interest
dat2 <- dplyr::select(dat2,roost_site,bat_age_class,bat_sex, collection_date,
                     species, sampleid, any_picorna,virus, type, family,
                     month, year, day_of_year, specific_roost_site)

head(dat2)
unique(dat2$roost_site)

#get into dat2e form
dat2$collection_date <- as.Date(dat2$collection_date,format = "%m/%d/%y")

#change picorna to numeric
dat2$any_picorna <- as.numeric(dat2$any_picorna)
names(dat2)[names(dat2)=="bat_species"] <- "species"

#and make sure it is only 1 of the same sample type per each dat2e
#dat2.list<- dlply(dat2, .(sampleid))

#summarize into prevalence by species and epiwk
dat2 <- ddply(dat2, .(sampleid, species, day_of_year,year,virus,family, bat_age_class, bat_sex, specific_roost_site), summarise, N=length(any_picorna), pos=sum(any_picorna))

#get negatives and prevalence
dat2$neg= dat2$N-dat2$pos
dat2$prevalence <- dat2$pos/dat2$N

#simplify= the name of "host_genus_species" 
names(dat2)[names(dat2)=="host_genus_species"] <- "species"

#make the dat2a into factors
dat2$species<-factor(dat2$species, levels=c("Eidolon dupreanum", "Pteropus rufus", "Rousettus madagascariensis"))
dat2$virus<-factor(dat2$virus, levels=c("aparavirus", "bat picornavirus", "cheravirus", "cripavirus", "felisavirus",
                                              "hepatovirus","kobuvirus","kunsagivirus","mischivirus","nepovirus","picorna-like virus",
                                              "sapelovirus","sapovirus","teschovirus","tetnovirus"))
dat2$family<-factor(dat2$family, levels=c("picornaviridae", "caliciviridae", "secoviridae","iflaviridae","unclassified"))
dat2$bat_sex<-factor(dat2$bat_sex, levels=c("Male", "Female"))
dat2$bat_age_class<-factor(dat2$bat_age_class, levels=c("J", "A"))
dat2$year<-factor(dat2$year, levels=c("2018","2019"))




#individual predictors of pos neg status exploration
all<-gam(pos~s(day_of_year, bs="cc", k=7)+s(year,bs="re"), data=dat2, family = binomial)
all
summary(all)
plot(all)
#not sig, no effect of year

species<-gam(pos~s(day_of_year,by=species, bs="cc", k=7)+s(year,bs="re"), data=dat2, family = binomial)
species
summary(species)
plot(species)
#RM and PR sig, no effect of year

virus<-gam(pos~s(day_of_year,by=virus, bs="cc", k=7)+s(year,bs="re"), data=dat2, family = binomial)
virus
summary(virus)
plot(virus)
#none sig, no effect of year

family<-gam(pos~s(day_of_year, by=family, bs="cc", k=7)+s(year,bs="re"), data=dat2, family = binomial)
family
summary(family)
plot(family)
#not sig, no effect of year

sex<-gam(pos~s(day_of_year, by=bat_sex, bs="cc", k=7)+s(year,bs="re"), data=dat2, family = binomial)
sex
summary(sex)
plot(sex)
#male sig, no effect of year

age<-gam(pos~s(day_of_year, by=bat_age_class,bs="cc", k=7)+s(year,bs="re"), data=dat2, family=binomial)
age
summary(age)
plot(age)
#adult sig, no effect of year






#Make nicer looking plots for "by bat"


#all positives
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(all, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(all, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

allpos<-ggplot(dat2)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Picornavirales prevalence over time")

allpos

#subset by species
ED<-subset(dat2, species=="Eidolon dupreanum")
RM<-subset(dat2, species=="Rousettus madagascariensis")
PR<-subset(dat2, species=="Pteropus rufus")

EDgam<-gam(pos~s(day_of_year,bs="cc", k=7)+s(year,bs="re"), data=ED, family = binomial)
RMgam<-gam(pos~s(day_of_year, bs="cc", k=7)+s(year,bs="re"), data=RM, family = binomial)
PRgam<-gam(pos~s(day_of_year, bs="cc", k=7)+s(year,bs="re"), data=PR, family = binomial)

#eidolon
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(EDgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(EDgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

edup<-ggplot(ED)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Eidolon dupreanum")

edup


#pteropus
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(PRgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(PRgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

pruf<-ggplot(PR)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Pteropus rufus")

pruf


#rousettus
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(RMgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(RMgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

rmad<-ggplot(RM)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
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
Adult<-subset(dat2, bat_age_class=="A")
Juv<-subset(dat2, bat_age_class=="J")
Adultgam<-gam(pos~s(day_of_year,bs="cc", k=7)+s(year,bs="re"), data=Adult, family=binomial)
Juvgam<-gam(pos~s(day_of_year, bs="cc", k=7)+s(year,bs="re"), data=Juv, family=binomial)

Fem<-subset(dat2, bat_sex=="Female")
Male<-subset(dat2, bat_sex=="Male")
Femgam<-gam(pos~s(day_of_year,bs="cc", k=7)+s(year,bs="re"), data=Fem, family=binomial)
Malegam<-gam(pos~s(day_of_year,bs="cc", k=7)+s(year,bs="re"), data=Male, family=binomial)


#females
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(Femgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Femgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

females<-ggplot(Fem)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Females")
females


#males
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(Malegam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Malegam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

males<-ggplot(Male)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Males")
males


#adults
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(Adultgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Adultgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

adults<-ggplot(Adult)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Adults")
adults


#juveniles
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(Juvgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(Juvgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

juveniles<-ggplot(Juv)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Juveniles")
juveniles


#plot next to each other for analysis
sexage<-plot_grid(females, males,adults, juveniles, nrow=2)
sexage



# #subset by virus
# apara<-subset(dat2, virus=="aparavirus")
# batpicorna<-subset(dat2, virus=="bat picornavirus")
# chera<-subset(dat2, virus=="cheravirus")
# cripa<-subset(dat2, virus=="cripavirus")
# felisa<-subset(dat2, virus=="felisavirus")
# hep<-subset(dat2, virus=="hepatovirus")
# kobu<-subset(dat2, virus=="kobuvirus")
# kun<-subset(dat2, virus=="kunsagivirus")
# mischi<-subset(dat2, virus=="mischivirus")
# nepo<-subset(dat2, virus=="nepovirus")
# picornalike<-subset(dat2, virus=="picorna-like virus")
# sapelo<-subset(dat2, virus=="sapelovirus")
# sapo<-subset(dat2, virus=="sapovirus")
# tescho<-subset(dat2, virus=="teschovirus")
# tetno<-subset(dat2, virus=="tetnovirus")
# 
# 
# apargam<-gam(pos~s(day_of_year,bs="cc")+s(year,bs="re"), data=apara, family=binomial)
# batpicornagam<-gam(pos~s(day_of_year, bs="cc",k=4)+s(year,bs="re"), data=batpicorna, family = binomial)
# cheragam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=chera, family = binomial)
# cripagam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=cripa, family = binomial)
# felisagam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=felisa, family = binomial)
# hepgam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=hep, family = binomial)
# kobugam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=kobu, family = binomial)
# kungam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=kun, family = binomial)
# mischigam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=mischi, family = binomial)
# nepogam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=nepo, family = binomial)
# picornalikegam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=picornalike, family = binomial)
# sapelogam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=sapelo, family = binomial)
# sapogam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=sapo, family = binomial)
# teschogam<-gam(pos~s(day_of_year, bs="cc", k=3)+s(year,bs="re"), data=tescho, family = binomial)
# tetnogam<-gam(pos~s(day_of_year, bs="cc")+s(year,bs="re"), data=tetno, family = binomial)
# 
# #aparavirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(aparagam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(aparagam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# aparavirus<-ggplot(apara)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Aparavirus")
# aparavirus
# 
# #bat picornavirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(batpicornagam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(batpicornagam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# batpicornavirus<-ggplot(batpicorna)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Bat picornavirus")
# batpicornavirus
# 
# #cheravirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(cheragam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(cheragam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# cheravirus<-ggplot(chera)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Cheravirus")
# cheravirus
# 
# #cripavirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(cripagam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(cripagam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# cripavirus<-ggplot(cripa)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Cripavirus")
# cripavirus
# 
# #felisavirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(felisagam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(felisagam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# felisavirus<-ggplot(felisa)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Felisavirus")
# felisavirus
# 
# #hepatovirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(hepgam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(hepgam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# hepatovirus<-ggplot(hep)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Hepatovirus")
# hepatovirus
# 
# #kobuvirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(kobugam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(kobugam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# kobuvirus<-ggplot(kobu)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Kobuvirus")
# kobuvirus
# 
# #kunsagivirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(kungam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(kungam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# kunsagivirus<-ggplot(kun)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Kunsagivirus")
# kunsagivirus
# 
# #mischivirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(mischigam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(mischigam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# mischivirus<-ggplot(mischi)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Mischivirus")
# mischivirus
# 
# #nepovirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(nepogam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(nepogam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# nepovirus<-ggplot(nepo)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Nepovirus")
# nepovirus
# 
# #picorna-like virus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(picornalikegam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(picornalikegam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# picornalikevirus<-ggplot(picornalike)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Picorna-like virus")
# picornalikevirus
# 
# #sapelovirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(sapelogam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(sapelogam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# sapelovirus<-ggplot(sapelo)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Sapelovirus")
# sapelovirus
# 
# #sapovirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(sapogam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(sapogam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# sapovirus<-ggplot(sapo)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Sapovirus")
# sapovirus
# 
# #teschovirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(teschogam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(teschogam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# teschovirus<-ggplot(tescho)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Teschovirus")
# teschovirus
# 
# #tetnovirus
# pred.df=cbind.data.frame(day_of_year=1:365)
# pred.df$predicted_count<-predict.gam(tetnogam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(tetnogam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# 
# tetnovirus<-ggplot(tetno)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lightsteelblue1", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="thistle1", alpha=0.05)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="olivedrab1", alpha=0.01)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle("Tetnovirus")
# tetnovirus
# 
# #plot next to each other for analysis
# viruses<-plot_grid(aparavirus, batpicornavirus, cheravirus, cripavirus,felisavirus, hepatovirus,
#                    kobuvirus, kunsagivirus, mischivirus, nepovirus, picornalikevirus,sapelovirus,
#                    sapovirus, teschovirus, tetnovirus, nrow=4)
# viruses


#subset by viral family
picorna<-subset(dat2, family=="picornaviridae")
cali<-subset(dat2, family=="caliciviridae")
seco<-subset(dat2, family=="secoviridae")
ifla<-subset(dat2, family=="iflaviridae")
un<-subset(dat2, family=="unclassified")

picornagam<-gam(pos~s(day_of_year, bs="cc", k=7)+s(year,bs="re"), data=picorna, family = binomial)
caligam<-gam(pos~s(day_of_year, bs="cc", k=5)+s(year,bs="re"), data=cali, family = binomial)
ungam<-gam(pos~s(day_of_year, bs="cc", k=4)+s(year,bs="re"), data=un, family = binomial)

#picornaviridae
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(picornagam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(picornagam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

picornaviridae<-ggplot(picorna)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Picornaviridae")
picornaviridae

#caliciviridae
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(caligam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(caligam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

caliciviridae<-ggplot(cali)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Caliciviridae")
caliciviridae

#unclassified
pred.df=cbind.data.frame(day_of_year=1:365, year=2018)
pred.df$predicted_count<-predict.gam(ungam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(ungam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

unclassified<-ggplot(un)+
  geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
  geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
  geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
  theme_classic()+ylab("Prevalence")+xlab("Day of year")+
  geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
  geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
  scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
  ggtitle("Unclassified family")
unclassified

#plot together
families<-plot_grid(picornaviridae,caliciviridae,unclassified, nrow=1)
families



