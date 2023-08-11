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


#####################################################################
# Set wd to data on this computer. Also ID homewd, assuming that 
# Mada-GIS is cloned to the same series of sub-folders
homewd = "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus" 
basewd = "/Users/gwenddolenkettenburg/Desktop/Mada-GIS"
#mapwd = paste0(basewd, "/", "Mada_GIS")
setwd(paste0(homewd, "/", "sampling_map/"))

#import madagascar shapfile
name<- paste0(basewd, "/", "MDG-3/MDG_adm3.shp")
otl_file <- paste(name, sep="") 
orotl_shp <- st_read(otl_file)
#View(orotl_shp)  # Open attribute table
class(orotl_shp)

dat <- read.csv(file = paste0(homewd,"/metadata/demo_picornavirales_fecesurine_meta_map_2018_2019.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

#remove double entries
bat.double = subset(dat, sampleid=="AMB122" |sampleid=="AMB123" |sampleid=="AMB124" |sampleid=="AMB126" |sampleid=="AMB127" |sampleid=="AMB128" 
                    |sampleid=="AMB129" |sampleid=="AMB130" |sampleid=="AMB131" |sampleid=="AMB132" |sampleid=="AMB133" |sampleid=="AMB134" |sampleid=="AMB135" 
                    |sampleid=="AMB138" |sampleid=="AMB139" |sampleid=="AMB140" |sampleid=="AMB141" |sampleid=="AMB142" |sampleid=="AMB143" |sampleid=="AMB144" 
                    |sampleid=="AMB145" |sampleid=="AMB146" |sampleid=="AMB147" |sampleid=="AMB149" |sampleid=="AMB150" |sampleid=="AMB151" |sampleid=="AMB152" 
                    |sampleid=="AMB153" |sampleid=="AMB154" |sampleid=="AMB155" |sampleid=="AMB159" |sampleid=="AMB160" |sampleid=="AMB162" |sampleid=="AMB165" 
                    |sampleid=="AMB170" |sampleid=="ANGB125" |sampleid=="ANGB127" |sampleid=="ANGB128" |sampleid=="ANGB129" |sampleid=="ANGB130" |sampleid=="ANGB131" 
                    |sampleid=="ANGB132" |sampleid=="ANGB133" |sampleid=="ANGB135" |sampleid=="ANGB137" |sampleid=="ANGB138" |sampleid=="ANGB141" |sampleid=="ANGB142" 
                    |sampleid=="ANGB155" |sampleid=="ANGB158" |sampleid=="ANGB161" |sampleid=="KEL148" |sampleid=="KEL149" |sampleid=="KEL150" |sampleid=="KEL151" 
                    |sampleid=="KEL152" |sampleid=="KEL153" |sampleid=="KEL154" |sampleid=="KEL155" |sampleid=="KEL156" |sampleid=="KEL157" |sampleid=="KEL159" 
                    |sampleid=="KEL160" |sampleid=="KEL161" |sampleid=="KEL162" |sampleid=="KEL163" |sampleid=="KEL164" |sampleid=="KEL165" |sampleid=="KEL166" 
                    |sampleid=="KEL167" |sampleid=="KEL168" |sampleid=="KEL169" |sampleid=="KEL170" |sampleid=="KEL171" |sampleid=="KEL172" |sampleid=="KEL173"
                    |sampleid=="KEL175" |sampleid=="KEL176" |sampleid=="KEL177" |sampleid=="KEL178" |sampleid=="KEL179" |sampleid=="KEL182" |sampleid=="KEL184" 
                    |sampleid=="KEL186" |sampleid=="KEL187" |sampleid=="KEL192" |sampleid=="KEL193" |sampleid=="KEL194" |sampleid=="KEL196" |sampleid=="KEL197" 
                    |sampleid=="KEL198" |sampleid=="KEL200" |sampleid=="KEL202" |sampleid=="KEL203" |sampleid=="KEL204" |sampleid=="KEL205" |sampleid=="KEL206" 
                    |sampleid=="KEL208" |sampleid=="KEL209" |sampleid=="KEL210" |sampleid=="KEL211" |sampleid=="KEL212" |sampleid=="KEL216" |sampleid=="KEL220" 
                    |sampleid=="KEL223" |sampleid=="KEL226" |sampleid=="KEL227" |sampleid=="KEL228" |sampleid=="KEL230" |sampleid=="KEL232" |sampleid=="KEL233" 
                    |sampleid=="KEL234" |sampleid=="KEL236" |sampleid=="KEL245" |sampleid=="KEL247" |sampleid=="KEL248" |sampleid=="KEL250" |sampleid=="KEL267" 
                    |sampleid=="KEL268" |sampleid=="KEL269" |sampleid=="KEL270" |sampleid=="KEL271" |sampleid=="KEL272" |sampleid=="KEL273" |sampleid=="KEL275" 
                    |sampleid=="MIZ135" |sampleid=="MIZ136" |sampleid=="MIZ137" |sampleid=="MIZ138" |sampleid=="MIZ139" |sampleid=="MIZ140" |sampleid=="MIZ141" 
                    |sampleid=="MIZ142" |sampleid=="MIZ143" |sampleid=="MIZ144" |sampleid=="MIZ145" |sampleid=="MIZ147" |sampleid=="MIZ148" |sampleid=="MIZ149" 
                    |sampleid=="MIZ150" |sampleid=="MIZ151" |sampleid=="MIZ152" |sampleid=="MIZ153" |sampleid=="MIZ154" |sampleid=="MIZ155" |sampleid=="MIZ156" 
                    |sampleid=="MIZ158" |sampleid=="MIZ159" |sampleid=="MIZ160" |sampleid=="MIZ161" |sampleid=="MIZ162" |sampleid=="MIZ163" |sampleid=="MIZ165" 
                    |sampleid=="MIZ167" |sampleid=="MIZ168" |sampleid=="MIZ169" |sampleid=="MIZ170" |sampleid=="MIZ171" |sampleid=="MIZ172" |sampleid=="MIZ173" 
                    |sampleid=="MIZ174" |sampleid=="MIZ175" |sampleid=="MIZ178" |sampleid=="MIZ179" |sampleid=="MIZ181" |sampleid=="MIZ184" |sampleid=="MIZ185" 
                    |sampleid=="MIZ186" |sampleid=="MIZ188" |sampleid=="MIZ190" |sampleid=="MIZ195" |sampleid=="MIZ207" |sampleid=="MIZ208" |sampleid=="MIZ275" 
                    |sampleid=="MIZ281" |sampleid=="MIZ298" |sampleid=="MIZ300" |sampleid=="MIZ305a" |sampleid=="MIZ305b" |sampleid=="MLB001" |sampleid=="MLB004" 
                    |sampleid=="MLB005" |sampleid=="MLB006")
unique(bat.double$any_picorna)


#and rank by rough age, by total positives
unique(dat$young_of_year)
bat.double$bat_age_class <- bat.double$bat_age_class
bat.double$bat_age_class[bat.double$bat_age_class=="P" | bat.double$bat_age_class=="L"] <- "A"
bat.double$bat_age_class[bat.double$bat_age_class=="NL" | bat.double$young_of_year=="no"] <- "A"
bat.double$bat_age_class[bat.double$young_of_year=="yes"] <- "J"

#by sampleID
unique(dat$young_of_year)
dat$bat_age_class <- dat$bat_age_class
dat$bat_age_class[dat$bat_age_class=="P" | dat$bat_age_class=="L"] <- "A"
dat$bat_age_class[dat$bat_age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$bat_age_class[dat$young_of_year=="yes"] <- "J"

#get into date form
bat.double$collection_date <- as.Date(bat.double$collection_date,format = "%m/%d/%y")

dat$collection_date <- as.Date(dat$collection_date,format = "%m/%d/%y")


#check sites
unique(dat$roost_site)

unique(bat.double$roost_site)



#get the date of the first day of every week
bat.double$epiwk <- cut(bat.double$collection_date, "week")
bat.double$epiwk <- as.Date(as.character(bat.double$epiwk))

dat$epiwk <- cut(dat$collection_date, "week")
dat$epiwk <- as.Date(as.character(dat$epiwk))

#change picorna to numeric
dat$any_picorna <- as.numeric(dat$any_picorna)

names(dat)[names(dat)=="bat_species"] <- "species"


bat.double$any_picorna <- as.numeric(bat.double$any_picorna)

names(dat)[names(dat)=="bat_species"] <- "species"

#and make sure it is only 1 of the same sample type per each date

dat.list <- dlply(dat, .(sampleid))

#summarize into prevalence by species and epiwk
dat.sum <- ddply(dat, .(species,epiwk), summarise, N=length(any_picorna), pos=sum(any_picorna))

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

#and plot
#here's a vector assigning colors to each species
colz = c("Eidolon dupreanum"="coral2", "Pteropus rufus" = "cornflowerblue", "Rousettus madagascariensis" = "darkgoldenrod1" )
#names(dat.sum)[names(dat.sum)=="bat_age_class"] <- "age"
# dat.sum$age[dat.sum$age=="A"] <- "adult"
# dat.sum$age[dat.sum$age=="J"] <- "juvenile"
# 
# shapez = c("juvenile" = 17, "adult" = 16)


#by points
p1 <- ggplot(data=dat.sum) + #here is the dataset
  geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci, color=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=epiwk, y= prevalence, color=species, size=N)) +#here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  geom_line(aes(x=epiwk, y= prevalence, color=species))+
  ylab("Picornavirus prevalence") + #change the name of the y-axis
  xlab("Sampling date") +
  scale_color_manual(values=colz) + #assign the colors manually using the vector above
  #facet_grid(age~.) +
  theme_bw() + #some style features
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
        strip.background = element_rect(fill="white"))
        #axis.title.x = element_blank())  #more style features
p1

