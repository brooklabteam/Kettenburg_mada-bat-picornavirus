
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
all<-gam(cbind(pos,neg)~s(day_of_year,by=bat_sex, bs="cc", k=7)+
           s(day_of_year,by=bat_age_class, bs="cc", k=7)+
           s(day_of_year,by=species, bs="cc", k=7)+
           s(year,bs="re"), data=dat2, family = binomial)
all
summary(all)
plot(all)
#Looks like R. mada is the only significant thing at p=0.0310

#Look at breakdown of demographic factors, first subset the data by species
edup<-subset(dat2, species=="Eidolon dupreanum")
rmad<-subset(dat2, species=="Rousettus madagascariensis")
pruf<-subset(dat2, species=="Pteropus rufus")

edupgam<-gam(cbind(pos,neg)~s(day_of_year,by=bat_sex, bs="cc", k=7)+
           s(day_of_year,by=bat_age_class, bs="cc", k=7)+
           s(year,bs="re"), data=edup, family = binomial)
edupgam
summary(edupgam)
plot(edupgam)
# none sig, but male looks like there's a pattern

rmadgam<-gam(cbind(pos,neg)~s(day_of_year,by=bat_sex, bs="cc", k=7)+
               s(day_of_year,by=bat_age_class, bs="cc", k=7)+
               s(year,bs="re"), data=rmad, family = binomial)
rmadgam
summary(rmadgam)
plot(rmadgam)
# none sig

prufgam<-gam(cbind(pos,neg)~s(day_of_year,by=bat_sex, bs="cc", k=7)+
               s(day_of_year,by=bat_age_class, bs="cc", k=7)+
               s(year,bs="re"), data=pruf, family = binomial)
prufgam
summary(prufgam)
plot(prufgam)
# none sig


#breakdown species by sex
edup_male<-subset(edup, bat_sex=="Male")
rmad_male<-subset(rmad, bat_sex=="Male")
pruf_male<-subset(pruf, bat_sex=="Male")
edup_female<-subset(edup, bat_sex=="Female")
rmad_female<-subset(rmad, bat_sex=="Female")
pruf_female<-subset(pruf, bat_sex=="Female")

#breakdown species by age
edup_a<-subset(edup, bat_age_class=="A")
rmad_a<-subset(rmad, bat_age_class=="A")
pruf_a<-subset(pruf, bat_age_class=="A")
edup_j<-subset(edup, bat_age_class=="J")
rmad_j<-subset(rmad, bat_age_class=="J")
pruf_j<-subset(pruf, bat_age_class=="J")

#combine sex and age into separate datasets by species
edup_male_a<-subset(edup_male, bat_age_class=="A")
rmad_male_a<-subset(rmad_male, bat_age_class=="A")
pruf_male_a<-subset(pruf_male, bat_age_class=="A")
edup_female_a<-subset(edup_female, bat_age_class=="A")
rmad_female_a<-subset(rmad_female, bat_age_class=="A")
pruf_female_a<-subset(pruf_female, bat_age_class=="A")
edup_male_j<-subset(edup_male, bat_age_class=="J")
rmad_male_j<-subset(rmad_male, bat_age_class=="J")
pruf_male_j<-subset(pruf_male, bat_age_class=="J")
edup_female_j<-subset(edup_female, bat_age_class=="J")
rmad_female_j<-subset(rmad_female, bat_age_class=="J")
pruf_female_j<-subset(pruf_female, bat_age_class=="J")

#make a dataset by species for viral families picornaviridae, unclassified, and caliciviridae
edup_cal<-subset(edup, family=="caliciviridae")
rmad_cal<-subset(rmad, family=="caliciviridae")
pruf_cal<-subset(pruf, family=="caliciviridae")
edup_pic<-subset(edup, family=="picornaviridae")
rmad_pic<-subset(rmad, family=="picornaviridae")
pruf_pic<-subset(pruf, family=="picornaviridae")
edup_un<-subset(edup, family=="unclassified")
rmad_un<-subset(rmad, family=="unclassified")
pruf_un<-subset(pruf, family=="unclassified")

#check out individual predictors of positivity now with the new datasets
##sub out the factors you're interested in "virus, family, etc" and the datasets
test<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
           s(year,bs="re"), data=rmad_male, family = binomial)
test
summary(test)
plot(test)

#in edup male dataset, by age neither A nor J is sig
#in edup female dataset, by age neither A nor J is sig
#in rmad male dataset, by age neither A nor J is sig
#in rmad female dataset, by age neither A nor J is sig
#in pruf male dataset, by age neither A nor J is sig
#in pruf female dataset, by age neither A nor J is sig
#adding in virus family and virus species is still all not sig
# in edup neither A nor J is sig, Male has pattern but is not sig. fem not sig
# in rmad Male is sig at 0.0044, female is not, neither A nor J is sig
# in pruf neither A nor J is sig but adult has pattern, neither A nor J is sig
# in pruf adult dataset, neiter male nor female is sig
# in pruf juvenile dataset, neiter male nor female is sig
# in edup a, male is sig at 0.0209, female is not
# in edup j, reduced knots to 4 and neither male nor female is sig
# in rmad a, male is sig at 0.0378 and female is not
# in rmad j, neither male nor female is sig



# Make plots of the following for main fig: 
#       edup adult data male (sig 0.02) and female, 
#       rmad adult data ,ale (sig 0.03) and female
# Make plots of the following for supplemental fig:
#       edup juvenile male and female, reduce knots to 4
#       rmad juvenile male and female
#       pruf adult male and female
#       pruf juvenile male and female


#break it down even more
edup_male_a<-subset(edup_male, bat_age_class=="A")
rmad_male_a<-subset(rmad_male, bat_age_class=="A")
pruf_male_a<-subset(pruf_male, bat_age_class=="A")
edup_female_a<-subset(edup_female, bat_age_class=="A")
rmad_female_a<-subset(rmad_female, bat_age_class=="A")
pruf_female_a<-subset(pruf_female, bat_age_class=="A")
edup_male_j<-subset(edup_male, bat_age_class=="J")
rmad_male_j<-subset(rmad_male, bat_age_class=="J")
pruf_male_j<-subset(pruf_male, bat_age_class=="J")
edup_female_j<-subset(edup_female, bat_age_class=="J")
rmad_female_j<-subset(rmad_female, bat_age_class=="J")
pruf_female_j<-subset(pruf_female, bat_age_class=="J")

edup_male_a_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
            s(year,bs="re"), data=edup_male_a, family = binomial)
edup_male_a_gam
summary(edup_male_a_gam)
plot(edup_male_a_gam)
# not sig

rmad_male_a_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
                       s(year,bs="re"), data=rmad_male_a, family = binomial)
rmad_male_a_gam
summary(rmad_male_a_gam)
plot(rmad_male_a_gam)
# not sig

pruf_male_a_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
                       s(year,bs="re"), data=pruf_male_a, family = binomial)
pruf_male_a_gam
summary(pruf_male_a_gam)
plot(pruf_male_a_gam)
# not sig

edup_female_a_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
                       s(year,bs="re"), data=edup_female_a, family = binomial)
edup_female_a_gam
summary(edup_female_a_gam)
plot(edup_female_a_gam)
# not sig

rmad_female_a_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
                       s(year,bs="re"), data=rmad_female_a, family = binomial)
rmad_female_a_gam
summary(rmad_female_a_gam)
plot(rmad_female_a_gam)
# not sig, but close

pruf_female_a_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
                       s(year,bs="re"), data=pruf_female_a, family = binomial)
pruf_female_a_gam
summary(pruf_female_a_gam)
plot(pruf_female_a_gam)
# not sig

# edup_male_j_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
#                        s(year,bs="re"), data=edup_male_j, family = binomial)
# edup_male_j_gam
# summary(edup_male_j_gam)
# plot(edup_male_j_gam)
# not enough data for this

rmad_male_j_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
                       s(year,bs="re"), data=rmad_male_j, family = binomial)
rmad_male_j_gam
summary(rmad_male_j_gam)
plot(rmad_male_j_gam)
# not sig

pruf_male_j_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
                       s(year,bs="re"), data=pruf_male_j, family = binomial)
pruf_male_j_gam
summary(pruf_male_j_gam)
plot(pruf_male_j_gam)
# not sig

edup_female_j_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=4)+
                         s(year,bs="re"), data=edup_female_j, family = binomial)
edup_female_j_gam
summary(edup_female_j_gam)
plot(edup_female_j_gam)
# not sig

rmad_female_j_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
                         s(year,bs="re"), data=rmad_female_j, family = binomial)
rmad_female_j_gam
summary(rmad_female_j_gam)
plot(rmad_female_j_gam)
# not sig, but close

# pruf_female_j_gam<-gam(cbind(pos,neg)~s(day_of_year, bs="cc", k=7)+
#                          s(year,bs="re"), data=pruf_female_j, family = binomial)
# pruf_female_j_gam
# summary(pruf_female_j_gam)
# plot(pruf_female_j_gam)
# not rnough data for this


##Nothing is really significant so I won't use these analysis in the paper!


#Make nicer looking plots for "by bat"

# #eidolon
# pred.df=cbind.data.frame(day_of_year=1:365, year=2019)
# pred.df$predicted_count<-predict.gam(EDgam, newdata=pred.df, type = "response", se.fit = T)$fit
# pred.df$predicted_count_se <- predict.gam(EDgam, newdata=pred.df,type = "response",se.fit = T)$se
# pred.df$predicted_count_lci <- pred.df$predicted_count -1.96*pred.df$predicted_count_se
# pred.df$predicted_count_uci <- pred.df$predicted_count +1.96*pred.df$predicted_count_se
# pred.df$predicted_count_lci[pred.df$predicted_count_lci <0] <- 0
# 
# title<-expression(paste(italic("Eidolon dupreanum, "), "p=0.714"))
# edup<-ggplot(ED)+
#   geom_rect(aes(xmin=152, xmax=244, ymin=-Inf, ymax=Inf),fill="lemonchiffon2", alpha=0.5)+
#   geom_rect(aes(xmin=217, xmax=314, ymin=-Inf, ymax=Inf),fill="#CCCCFF", alpha=0.02)+
#   geom_rect(aes(xmin=314, xmax=365, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
#   geom_rect(aes(xmin=1, xmax=31, ymin=-Inf, ymax=Inf),fill="paleturquoise1", alpha=0.1)+
#   theme_classic()+ylab("Prevalence")+xlab("Day of year")+
#   geom_line(data=pred.df, aes(x=day_of_year, y=predicted_count), size=1) +
#   geom_ribbon(data=pred.df, aes(x=day_of_year, ymin=predicted_count_lci, ymax=predicted_count_uci), alpha=.3)+
#   theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),axis.title.y = element_text(size=10),axis.title.x = element_text(size=10))+
#   scale_x_continuous(breaks=c(0,100,200,300), labels=c("0", "100", "200", "300"))+
#   ggtitle(title)
# 
# edup
# 
# 
# 
# 
# #plot next to each other for analysis
# Fig2<-plot_grid(females, males,adults, juveniles, nrow=2, labels="AUTO", label_size = 23)
# Fig2

#export 10x8inch landscape PDF for Fig 2

