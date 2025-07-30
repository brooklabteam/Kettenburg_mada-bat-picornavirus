rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)
library(phylotools)
library(phytools)
library(phylobase)
library(cowplot)
library(ggtreeExtra)
library(viridis)
library(ggplotify)
library(patchwork)

###packages loaded
##########################################################################################################
##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_mada-bat-picornavirus/"
setwd(paste0(homewd,"/tanglegram"))

#Bat picorna tanglegram

tr1 <- read.tree("batpicorna_hosts_short.newick")
tr2 <- read.tree("batpicorna_fullname_align_fasttree_short.newick")

# assoc<-cbind(c("Rousettus_aegyptiacus",
#                "Rousettus_aegyptiacus",
#                "Rousettus_aegyptiacus",
#                "Rhinolophus_sinicus",
#                "Taphozous_melanopogon",
#                "Rousettus_madagascariensis",
#                "Rousettus_madagascariensis",
#                "Rousettus_madagascariensis",
#                "Rousettus_madagascariensis",
#                "Rousettus_madagascariensis",
#                "Rousettus_aegyptiacus",
#                "Rousettus_aegyptiacus",
#                "Rousettus_aegyptiacus"),
#              
#              c("PP711928.1_Rousettus_bat_picornavirus_Uganda" ,        
#                "PP711909.1_Rousettus_bat_picornavirus_Uganda"  ,       
#                 "PP711930.1_Rousettus_bat_picornavirus_Kenya" ,         
#                 "PP746000.1_Batpicornavirus_BtSY4_China" ,              
#                 "PP745912.1_Bat_picornavirus_7_China"  ,                
#                 "OQ818325_R_madagascariensis_picornavirus_1_Madagascar",
#                 "OQ818328_R_madagascariensis_picornavirus_1_Madagascar",
#                 "PP766472_R_madagascariensis_picornavirus_3_Madagascar",
#                 "PP766469_R_madagascariensis_picornavirus_3_Madagascar",
#                 "PV788825_R_madagascariensis_picornavirus_4_Madagascar",
#                 "PP711945.1_Rousettus_bat_picornavirus_Uganda" ,        
#                 "PP711912.1_Rousettus_bat_picornavirus_Uganda",         
#                 "PP711913.1_Rousettus_bat_picornavirus_Kenya"))

assoc<-cbind(c("R_aegyptiacus",
               "R_aegyptiacus",
               "R_aegyptiacus",
               "R_sinicus",
               "T_melanopogon",
               "R_madagascariensis",
               "R_madagascariensis",
               "R_madagascariensis",
               "R_madagascariensis",
               "R_madagascariensis",
               "R_aegyptiacus",
               "R_aegyptiacus",
               "R_aegyptiacus"),
             
             c("PP711928.1" ,        
               "PP711909.1"  ,       
               "PP711930.1" ,         
               "PP746000.1" ,              
               "PP745912.1"  ,                
               "OQ818325*",
               "OQ818328*",
               "PP766472*",
               "PP766469*",
               "PV788825*",
               "PP711945.1" ,        
               "PP711912.1",         
               "PP711913.1"))


batpicorna<-cophylo(tr1,tr2,assoc,rotate=TRUE)
plot(batpicorna, fsize=0.6, gap=2)


#Bat sapelo tanglegram

tr1 <- read.tree("sapelo_hosts_short.newick")
tr2 <- read.tree("batsapelo_fullname_align_fasttree_short.newick")

# assoc<-cbind(c("Eidolon_helvum",
             #   "Eidolon_dupreanum",
             #   "Eidolon_helvum",
             #   "Eidolon_helvum",
             #   "Eidolon_dupreanum",
             #   "Rousettus_aegyptiacus",
             #   "Rousettus_madagascariensis",
             #   "Eonycteris_spelaea",
             #   "Eonycteris_spelaea",
             #   "Eonycteris_spelaea",
             #   "Rousettus_leschenaultii"),
             # 
             # c("PP711921.1_Eidolon_bat_sapelovirus_Kenya" ,           
             #   "OQ818321_E_dupreanum_sapelovirus_2_Madagascar"    ,   
             #    "PP711943.1_Eidolon_bat_sapelovirus_Kenya"  ,          
             #    "NC_033820.1_Batsapelovirus_Cameroon"  ,               
             #   "OQ818320_E_dupreanum_sapelovirus_1_Madagascar" ,      
             #    "PP711911.1_Rousettus_bat_sapelovirus_Kenya" ,         
             #    "OQ818329_R_madagascariensis_sapelovirus_1_Madagascar",
             #   "OR951326.1_Pteropodidae_bat_sapelovirus_China" ,      
             #  "OR951327.1_Pteropodidae_bat_sapelovirus_China"  ,     
             #    "OR951325.1_Pteropodidae_bat_sapelovirus_China"  ,     
             #    "OR951332.1_Pteropodidae_bat_sapelovirus_China" ))
assoc<-cbind(c("E_helvum",
               "E_dupreanum",
               "E_helvum",
               "E_helvum",
               "E_dupreanum",
               "R_aegyptiacus",
               "R_madagascariensis",
               "E_spelaea",
               "E_spelaea",
               "E_spelaea",
               "R_leschenaultii"),
             
             c("PP711921.1" ,           
               "OQ818321*"    ,   
               "PP711943.1"  ,          
               "NC_033820.1"  ,               
               "OQ818320*" ,      
               "PP711911.1" ,         
               "OQ818329*",
               "OR951326.1" ,      
               "OR951327.1"  ,     
               "OR951325.1"  ,     
               "OR951332.1" ))

sapelo<-cophylo(tr1,tr2,assoc,rotate=TRUE)
plot(sapelo,fsize=0.55)



#Bat sapo tanglegram

tr1 <- read.tree("sapo_hosts_short.newick")
tr2 <- read.tree("batsapo_fullname_align_fasttree_short.newick")

# assoc<-cbind(c("Cynopterus_sphinx",
#                "Rousettus_leschenaultii" ,
#                "Rousettus_aegyptiacus",
#                "Rousettus_madagascariensis",
#                "Rousettus_leschenaultii" ,
#                "Rousettus_aegyptiacus",
#                "Rousettus_madagascariensis",
#                "Rousettus_leschenaultii",
#                "Rousettus_madagascariensis",
#                "Rousettus_aegyptiacus",
#                "Rousettus_aegyptiacus",
#                "Rousettus_madagascariensis",
#                "Rousettus_madagascariensis",
#                "Rousettus_madagascariensis",
#                "Eidolon_helvum",
#                "Eidolon_helvum",
#                "Eonycteris_spelaea",
#                "Eidolon_helvum",
#                "Eidolon_helvum",
#                "Eidolon_helvum",
#                "Eidolon_helvum",
#                "Eidolon_helvum"),
#              
#              c(  "OP963622.1_Bat_sapovirus_China",                    
#                  "OQ709197.1_Bat_sapovirus_China"   ,                 
#                  "PP712008.1_Rousettus_bat_calicivirus_Kenya"   ,     
#                  "PV788824_R_madagascariensis_sapovirus_4_Madagascar",
#                "OP963623.1_Bat_sapovirus_China"     ,               
#                  "PP712015.1_Rousettus_bat_calicivirus_Kenya"        ,
#                  "PP766470_R_madagascariensis_sapovirus_2_Madagascar",
#                  "OR951145.1_Bat_sapovirus_China"    ,                
#                  "OQ818347_R_madagascariensis_sapovirus_2_Madagascar",
#                  "PP712001.1_Rousettus_bat_calicivirus_Kenya"        ,
#                 "PP712004.1_Rousettus_bat_calicivirus_Kenya"        ,
#                  "OQ818340_E_dupreanum_sapovirus_2_Madagascar"       ,
#                  "OQ818319_E_dupreanum_sapovirus_1_Madagascar"       ,
#                  "PP766459_E_dupreanum_sapovirus_1_Madagascar"       ,
#                  "KX759619.1_Batsapovirus_Cameroon"                  ,
#                  "KX759623.1_Batsapovirus_Cameroon"                  ,
#                 "OR951138.1_Batsapovirus_China"                     ,
#                 "KX759620.1_Batsapovirus_Cameroon"                  ,
#                 "NC_033776.1_Batsapovirus_Cameroon"                 ,
#                  "KX759622.1_Batsapovirus_Cameroon"                  ,
#                  "KX759618.1_Batsapovirus_Cameroon"                  ,
#                  "KX759621.1_Batsapovirus_Cameroon"))

assoc<-cbind(c("C_sphinx",
               "R_leschenaultii" ,
               "R_aegyptiacus",
               "R_madagascariensis",
               "R_leschenaultii" ,
               "R_aegyptiacus",
               "R_madagascariensis",
               "R_leschenaultii",
               "R_madagascariensis",
               "R_aegyptiacus",
               "R_aegyptiacus",
               "E_dupreanum",
               "E_dupreanum",
               "E_dupreanum",
               "E_helvum",
               "E_helvum",
               "E_spelaea",
               "E_helvum",
               "E_helvum",
               "E_helvum",
               "E_helvum",
               "E_helvum"),
             
             c(  "OP963622.1",                    
                 "OQ709197.1"   ,                 
                 "PP712008.1"   ,     
                 "PV788824*",
                 "OP963623.1"     ,               
                 "PP712015.1"        ,
                 "PP766470*",
                 "OR951145.1"    ,                
                 "OQ818347*",
                 "PP712001.1"        ,
                 "PP712004.1"        ,
                 "OQ818340*"       ,
                 "OQ818319*"       ,
                 "PP766459*"       ,
                 "KX759619.1"                  ,
                 "KX759623.1"                  ,
                 "OR951138.1"                     ,
                 "KX759620.1"                  ,
                 "NC_033776.1"                 ,
                 "KX759622.1"                  ,
                 "KX759618.1"                  ,
                 "KX759621.1"))

sapo<-cophylo(tr1,tr2,assoc,rotate=TRUE)
plot(sapo,fsize=0.55)

#Bat tescho tanglegram

tr1 <- read.tree("tescho_hosts_short.newick")
tr2 <- read.tree("battescho_fullname_align_fasttree_short.newick")

# assoc<-cbind(c("Rousettus_leschenaultii" ,
#                "Rousettus_madagascariensis",
#                "Rousettus_madagascariensis",
#                "Rousettus_aegyptiacus",
#                "Rousettus_madagascariensis",
#                "Eonycteris_spelaea",
#                "Rousettus_leschenaultii" ,
#                "Rousettus_leschenaultii" ,
#                "Rousettus_aegyptiacus",
#                "Eidolon_helvum",
#                "Eidolon_dupreanum"),
#              
#              c("OR951335.1_Pteropodidae_bat_teschovirus_China",       
#                "OQ818324_R_madagascariensis_teschovirus_2_Madagascar",
#                "PV788826_R_madagascariensis_teschovirus_2_Madagascar",
#                 "PP711934.1_Rousettus_bat_teschovirus_Kenya",          
#                 "OQ818323_R_madagascariensis_teschovirus_1_Madagascar",
#                 "OR951324.1_Pteropodidae_bat_teschovirus_China" ,      
#                 "OR951333.1_Pteropodidae_bat_teschovirus_China",       
#                 "OR951334.1_Pteropodidae_bat_teschovirus_China",       
#                 "PP711948.1_Rousettus_bat_teschovirus_Uganda",         
#                 "KX420938.1_Teschovirus_sp_Saudi_Arabia",              
#                 "OQ818318_E_dupreanum_teschovirus_1_Madagascar" ))

assoc<-cbind(c("R_leschenaultii" ,
               "R_madagascariensis",
               "R_madagascariensis",
               "R_aegyptiacus",
               "R_madagascariensis",
               "E_spelaea",
               "R_leschenaultii" ,
               "R_leschenaultii" ,
               "R_aegyptiacus",
               "E_helvum",
               "E_dupreanum"),
             
             c("OR951335.1",       
               "OQ818324*",
               "PV788826*",
               "PP711934.1",          
               "OQ818323*",
               "OR951324.1" ,      
               "OR951333.1",       
               "OR951334.1",       
               "PP711948.1",         
               "KX420938.1",              
               "OQ818318*" ))

tescho<-cophylo(tr1,tr2,assoc,rotate=TRUE)
plot(tescho,fsize=0.55)


#hepato tanglegram

tr1 <- read.tree("hepato_hosts_short.newick")
tr2 <- read.tree("hepato_fullname_align_fasttree_short.newick")

# assoc<-cbind(c("Sorex_araneus" ,
#                "Miniopterus_manavi",
#                "Eptesicus_fuscus",
#                "Hipposideros_armiger",
#                "Hipposideros_larvatus",
#                "Artibeus_planirostris",
#                "Tupaia_belangeri",
#                "Eidolon_helvum",
#                "Eidolon_dupreanum",
#                "Eidolon_dupreanum",
#                "Eidolon_dupreanum"),
#              
#              c( "KT452661.1_shrew_hepatovirus_Germany",         
#                  "NC_038313.1_bat_hepatovirus_Madagascar"  ,     
#                 "OM302498.1_e_fuscus_USA"   ,                   
#                  "MG559674.1_hepatovirus_sp_China"   ,           
#                  "OR951271.1_hipposideros_bat_hepatovirus_China",
#                 "OR367421.1_hepatovirus_sp_Brazil"  ,           
#                 "NC_028981.1_Tupaia_hepatovirus_A_China"   ,    
#                 "NC_028366.1_Hepatovirus_H2_Ghana"   ,          
#                 "PP766455_E_dupreanum_hepatovirus_Madagascar"  ,
#                  "PP766457_E_dupreanum_hepatovirus_Madagascar" , 
#                  "OQ818337_E_dupreanum_hepatovirus_Madagascar"   ))

assoc<-cbind(c("S_araneus" ,
               "M_manavi",
               "E_fuscus",
               "H_armiger",
               "H_larvatus",
               "A_planirostris",
               "T_belangeri",
               "E_helvum",
               "E_dupreanum",
               "E_dupreanum",
               "E_dupreanum"),
             
             c( "KT452661.1",         
                "NC_038313.1"  ,     
                "OM302498.1"   ,                   
                "MG559674.1"   ,           
                "OR951271.1",
                "OR367421.1"  ,           
                "NC_028981.1"   ,    
                "NC_028366.1"   ,          
                "PP766455*"  ,
                "PP766457*" , 
                "OQ818337*"   ))

hepato<-cophylo(tr1,tr2,assoc,rotate=TRUE)
plot(hepato,fsize=0.55)


#mischi tanglegram

tr1 <- read.tree("mischi_hosts_short.newick")
tr2 <- read.tree("mischi_fullname_align_fasttree_short.newick")

# assoc<-cbind(c("Canis_lupus",
#                "Hipposideros_gigas",
#                "Pteropus_rufus",
#                "Miniopterus_pusillus",
#                "Miniopterus_schreibersii",
#                "Miniopterus_pusillus",
#                "Miniopterus_schreibersii",
#                "Miniopterus_schreibersii",
#                "Miniopterus_schreibersii",
#                "Miniopterus_schreibersii",
#                "Miniopterus_schreibersii",
#                "Miniopterus_schreibersii"),
#              
#              c( "NC_075428.1_Mischivirus_D_USA",          
#                "NC_026470.1_Mischivirus_C_DRC" ,         
#                 "OQ818316_P_rufus_mischivirus_Madagascar",
#                 "OR951306.1_Bat_mischivirus_5_China"  ,   
#                 "NC_043072.1_Mischivirus_B_Algeria"  ,    
#                 "OR951299.1_Mischivirus_A_China" ,        
#                 "NC_034381.1_Mischivirus_A_China" ,       
#                 "OR867092.1_Mischivirus_sp_China"   ,     
#                 "OR951342.1_Bat_mischivirus_4_China" ,    
#                 "OR951362.1_Bat_mischivirus_4_China" ,    
#                "OR951339.1_Bat_mischivirus_4_China",     
#                 "OR951347.1_Bat_mischivirus_4_China" ))

assoc<-cbind(c("C_lupus",
               "H_gigas",
               "P_rufus",
               "M_pusillus",
               "M_schreibersii",
               "M_pusillus",
               "M_schreibersii",
               "M_schreibersii",
               "M_schreibersii",
               "M_schreibersii",
               "M_schreibersii",
               "M_schreibersii"),
             
             c( "NC_075428.1",          
                "NC_026470.1" ,         
                "OQ818316*",
                "OR951306.1"  ,   
                "NC_043072.1"  ,    
                "OR951299.1" ,        
                "NC_034381.1" ,       
                "OR867092.1"   ,     
                "OR951342.1" ,    
                "OR951362.1" ,    
                "OR951339.1",     
                "OR951347.1" ))

mischi<-cophylo(tr1,tr2,assoc,rotate=TRUE)
plot(mischi,fsize=0.55)










