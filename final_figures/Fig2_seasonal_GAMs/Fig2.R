
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
all<-gam(pos~s(day_of_year,by=bat_sex, bs="cc", k=7)+
           s(day_of_year,by=bat_age_class, bs="cc", k=7)+
           s(day_of_year,by=species, bs="cc", k=7)+
           s(year,bs="re"), data=dat2, family = binomial)
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
  ggtitle("Picornavirales prevalence over time, p=0.113")

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

title<-expression(paste(italic("Eidolon dupreanum, "), "p=0.714"))
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
  ggtitle(title)

edup


#pteropus
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(PRgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(PRgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

title<-expression(paste(italic("Pteropus rufus, "), "p=0.034*"))
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
  ggtitle(title)

pruf


#rousettus
pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
pred.df$predicted_count<-predict.gam(RMgam, newdata=pred.df, type = "response", se.fit = T)$fit
pred.df$predicted_count_se <- predict.gam(RMgam, newdata=pred.df,type = "response",se.fit = T)$se
pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0

title<-expression(paste(italic("Rousettus madagascariensis, "), "p=0.006**"))
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
  ggtitle(title)
rmad

#all three next to each other for analysis
species<-plot_grid(edup, pruf, rmad, nrow=1, labels="AUTO", label_size = 23)
species

#export 15x5 inch PDF landscape as a supplemental figure

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
  ggtitle("Females, p=0.782")
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
  ggtitle("Males, p=0.026*")
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
  ggtitle("Adults, p=0.090")
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
  ggtitle("Juveniles, p=0.761")
juveniles


#plot next to each other for analysis
Fig2<-plot_grid(females, males,adults, juveniles, nrow=2, labels="AUTO", label_size = 23)
Fig2

#export 10x8inch landscape PDF for Fig 2

