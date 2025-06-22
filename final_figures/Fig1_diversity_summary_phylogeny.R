#This script is for figure 1, which is a master phylogeny using all the viral genera identified for the 
#RDRP gene to create an infographic, then additional genus specific phylogenies to go with it. It also has a panel
#to show the diversity of viral species in each sampling session

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
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_mada-bat-picornavirus/"

##Add a plot showing the shared number of viruses in population during each sampling session
dat <- read.csv(file = paste0(homewd,"/metadata/demo_data_indiv_pos_heatmap_genus.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

#pick order for the labels
dat$genus <- factor(dat$genus, levels = c("Cardio.","Hepato.","Kobu.","Kunsagi.","Mischi.", "Bat picorna.",
                                          "Sapelo.","Tescho.", "Sapo."))   
dat$virus <- factor(dat$virus, levels = c("Eidolon dupreanum cardiovirus", "Eidolon dupreanum hepatovirus", "Eidolon dupreanum kobuvirus", "Eidolon dupreanum kobuvirus 2",
                                          "Eidolon dupreanum kunsagivirus", "Pteropus rufus mischivirus", "Rousettus madagascariensis picornavirus 1", "Rousettus madagascariensis picornavirus 2",
                                          "Rousettus madagascariensis picornavirus 3","Rousettus madagascariensis picornavirus 4", "Eidolon dupreanum sapelovirus 1","Eidolon dupreanum sapelovirus 2", "Rousettus madagascariensis sapelovirus 1",
                                          "Eidolon dupreanum teschovirus 1","Rousettus madagascariensis teschovirus 1","Rousettus madagascariensis teschovirus 2","Eidolon dupreanum sapovirus 1",
                                          "Eidolon dupreanum sapovirus 2", "Eidolon dupreanum sapovirus 3", "Eidolon dupreanum sapovirus 4", "Rousettus madagascariensis sapovirus 1",
                                          "Rousettus madagascariensis sapovirus 2", "Rousettus madagascariensis sapovirus 3", "Rousettus madagascariensis sapovirus 4"))   

#pick colors for virus genera
# genuscolz<- c("Cardio."="#F8766D","Hepato."="#D89000","Kobu."="#A3A500","Kunsagi."="#39B600","Mischi."="#00BF7D",
#               "Sapelo."="#00BFC4","Sapo."="#00B0F6","Tescho."="#E76BF3","Bat picorna."="#9590FF",
#               "Alpha."="black")
genuscolz<- c("Cardiovirus"="#F8766D","Hepatovirus"="#D89000","Kobuvirus"="#A3A500","Kunsagivirus"="#39B600","Mischivirus"="#00BF7D",
              "Sapelovirus"="#00BFC4","Sapovirus"="#00B0F6","Teschovirus"="#E76BF3","Bat picornavirus"="#9590FF",
              "Alphavirus"="black")

library(scales)
hex_codes2 <- hue_pal()(30)      
show_col(hex_codes2)

viruscolz<-c("Eidolon dupreanum cardiovirus"="#F8766D", "Eidolon dupreanum hepatovirus"="#D89000", "Eidolon dupreanum kobuvirus"="#A3A500", "Eidolon dupreanum kobuvirus 2"="#A3A560",
             "Eidolon dupreanum kunsagivirus"="#39B600", "Pteropus rufus mischivirus"="#00BF7D", "Rousettus madagascariensis picornavirus 1"="#9590FF", "Rousettus madagascariensis picornavirus 2"="#B983FF",
             "Rousettus madagascariensis picornavirus 3"="mediumpurple4", "Rousettus madagascariensis picornavirus 4"="thistle4","Eidolon dupreanum sapelovirus 1"="aquamarine1","Eidolon dupreanum sapelovirus 2"="aquamarine3", "Rousettus madagascariensis sapelovirus 1"="aquamarine4",
             "Eidolon dupreanum teschovirus 1"="palevioletred1","Rousettus madagascariensis teschovirus 1"="palevioletred3","Rousettus madagascariensis teschovirus 2"="palevioletred4","Eidolon dupreanum sapovirus 1"="lightskyblue1",
             "Eidolon dupreanum sapovirus 3"="lightskyblue3", "Eidolon dupreanum sapovirus 4"="lightskyblue4", "Rousettus madagascariensis sapovirus 1"="#619CFF",
             "Rousettus madagascariensis sapovirus 2"="royalblue3", "Rousettus madagascariensis sapovirus 3"="royalblue4", "Rousettus madagascariensis sapovirus 4"="navy")


dat$bat_species<-factor(dat$bat_species, levels=c("Pteropus rufus","Eidolon dupreanum", "Rousettus madagascariensis"))
dat$roost_site<-factor(dat$roost_site, levels=c("Angavokely/Angavobe","Ambakoana","Maromizaha"))
dat$sampling_session<-factor(dat$sampling_session)

#Subset because Pteropus only has one sample pos with only one virus
dat<-subset(dat, bat_species!="Pteropus rufus")

library(ggh4x)

#virus and genus together
# p4<-ggplot(dat) +
#   geom_bar(aes(x = sampling_date, y = num_virus, fill = virus),
#            position = "stack",
#            stat = "identity") +
#   #facet_nested(roost_site+genus~., scales="free", space="free")+
#   facet_nested(genus~roost_site, scales="free", space="free",
#                 nest_line = element_line(color="white"), solo_line = TRUE)+
#   scale_fill_manual(values=viruscolz)+
#   
#   labs(
#     x = "Sampling date",
#     y= "Number of sequences",
#     fill="Virus species",
#     title="")+
#   theme_linedraw()+
#   scale_y_continuous(n.breaks = 2)+
#   guides(fill = guide_legend(ncol = 1))+
#   theme(plot.margin = margin(0, 1, 0, 30, "pt"),
#         plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.y = element_text(size=10),
#         #axis.text.x = element_text(size=10),
#         axis.text.x = element_text(angle = 90, size=10),
#         legend.text = element_text(size=8, face="italic"),
#         legend.title = element_text(size=9),
#         legend.position = "right")
# 
# p4


sampling<-ggplot(dat,aes(x = sampling_date, y = num_genus, fill = genus_full, label=num_virus)) +
  geom_bar(
    position = "stack",
    stat = "identity") +
  #geom_text(size = 3, position = position_stack(vjust = 0.5))+
  facet_nested(.~roost_site, scales="free", space="free",
               nest_line = element_line(color="white"), solo_line = TRUE)+
  scale_fill_manual(values=genuscolz)+
  
  labs(
    x = "Sampling date",
    y= "Number of unique viral species per genus",
    fill="Virus genus",
    title="")+
  theme_linedraw()+
  scale_y_continuous(n.breaks = 4)+
  guides(fill = guide_legend(ncol = 1))+
  theme(plot.margin = margin(0, 1, 0, 30, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        axis.text.y = element_text(size=10),
        #axis.text.x = element_text(size=10),
        axis.text.x = element_text(angle = 90, size=10),
        legend.text = element_text(size=8, face="italic"),
        legend.title = element_text(size=9),
        legend.position = "right")

sampling

#all tree code
##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_mada-bat-picornavirus/"
setwd(paste0(homewd,"/IQtree_phylogenies/master_phylo_fig"))

#load the tree and root it
tree <-  read.tree("polymerase_summary_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#load tree data prepared from elsewhere
dat <- read.csv(("polymerase_summary_metadata_circ.csv"), header = T, stringsAsFactors = F)

head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #296
length(tree$tip.label) #296

#check subgroup names
unique(dat$Genus)

#pick order for the labels
dat$Genus <- factor(dat$Genus, levels = c("Cardiovirus","Hepatovirus","Kobuvirus","Kunsagivirus","Mischivirus",
                                          "Sapelovirus","Sapovirus","Teschovirus","Unclassified bat picornavirus",
                                          "Alphavirus"))   
#pick colors for virus genera
genuscolz<- c("Cardiovirus"="#F8766D","Hepatovirus"="#D89000","Kobuvirus"="#A3A500","Kunsagivirus"="#39B600","Mischivirus"="#00BF7F",
              "Sapelovirus"="#00BFC4","Sapovirus"="#00B0F6","Teschovirus"="#E76BF3","Unclassified bat picornavirus"="#9590FF",
              "Alphavirus"="black")

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  scale_color_manual(values=genuscolz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now check that you don't have NAs and blanks throughout the component parts of the name
#you also can edit the original csv file to replace these blanks if you can find the correct info

dat$Isolate  #some are blank, so convert those to NA
dat$Isolate[dat$Isolate==""] <- NA
dat$Accession #all good

dat$source
dat$source[dat$source==""] <- NA
dat$Country #some are blank, so convert those to NA
dat$Country[dat$Country==""] <- NA
dat$Collection_Date #these are messy, some are years and some are full dates. I just want years, will manually fix


#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Genus[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Genus[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[ !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Genus[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


#look at dat$new_label
dat$new_label #they all look great

#make sure to sort in order
colnames(dat)[colnames(dat)=="Accession"]="old_tip_label"
tree.dat <- data.frame(old_tip_label=rooted.tree$tip.label, num =1:length(rooted.tree$tip.label))
head(tree.dat)
head(dat)
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = F, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label

#making separate datasets for heatmaps to go next to tree
extra<-data.frame(id=tree.dat$tip_label,Novel_contigs=tree.dat$num_genus_novel_contigs,
                  Novel_reads_log10=tree.dat$num_genus_novel_reads_log,Region=tree.dat$region,
                  Host_class=tree.dat$host_class)
#rownames(extra) <- rooted.tree$tip.label

contig<-data.frame(id=tree.dat$tip_label,Novel_contigs=tree.dat$num_genus_novel_contigs)
#rownames(contig) <- rooted.tree$tip.label

reads<-data.frame(id=tree.dat$tip_label,Novel_reads_log10=tree.dat$num_genus_novel_reads_log)
#rownames(reads) <- rooted.tree$tip.label

region<-data.frame(id=tree.dat$tip_label,Region=tree.dat$region)
#rownames(region) <- rooted.tree$tip.label

host<-data.frame(id=tree.dat$tip_label,Host_class=tree.dat$host_class)
#rownames(host) <- rooted.tree$tip.label

novel<-data.frame(id=tree.dat$tip_label,Seq_type=tree.dat$Seq_type)
#rownames(novel) <- rooted.tree$tip.label

#make real tree.dat file for putting the tip points on the tree
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Seq_type, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new

cbind(tree.dat$old_tip_label, rooted.tree$tip.label)

#check out the labels
tree.dat$tip_label#all look good

#assign some stuff for shapes
tree.dat$Host[tree.dat$Host==0] <- "Non-bat host"
tree.dat$Host[tree.dat$Host==1] <- "Bat host"
tree.dat$Host <- as.factor(tree.dat$Host)
tree.dat$Seq_type[tree.dat$Seq_type==0] <- "Reference seq"
tree.dat$Seq_type[tree.dat$Seq_type==1] <- "Novel seq"
tree.dat$Seq_type<-as.factor(tree.dat$Seq_type)
novel$Seq_type[novel$Seq_type==0] <- "Reference seq"
novel$Seq_type[novel$Seq_type==1] <- "Novel seq"
#novel$Seq_type<-as.factor(novel$Seq_type)

shapez = c("Bat host" =  17, "Non-bat host" = 19, "Reference seq"=19,"Novel seq"=17)
colz2 = c('1' =  "tomato1", '0' = "white")

circ<-ggtree(rooted.tree, layout="circular")

##Get the clade numbers so we can label
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)

##base tree
p1 <- ggtree(rooted.tree, layout="fan", size=0.5) %<+% tree.dat +
  geom_tippoint(aes(color=Genus, shape=Host), size=2,stroke=0,show.legend = T) +
  #scale_color_manual(values=genuscolz)+
  scale_shape_manual(values=shapez) +
  guides(colour = "none", shape = guide_legend(order = 1))+
  theme(#legend.position = c(0.5,0.59), #keep this one in case we want the legend within the plot
        #legend.position = c(0.97,0.59), #right side
        legend.position = c(0.52,0),
        legend.margin = margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.25, "cm"),
        legend.direction = "vertical",
        legend.box = "horizontal") +
  xlim(c(0,9)) 
p1

library(scales)
hex_codes2 <- hue_pal()(10)      
show_col(hex_codes2)

#rotate tree a bit
p1<-rotate_tree(p1, 30)
p1

#add clade labels
p1.1 <- p1 +
  geom_cladelabel(node = 403, label = 'italic(Kobuvirus)', parse=TRUE,offset=2.5, offset.text=0.4, fontsize=4, angle=310, hjust=0.6,align = TRUE, color="#A3A500") +
  geom_cladelabel(node = 461, label = "italic(Kunsagivirus)", parse=TRUE, offset=2.5, offset.text=0.4,  fontsize=4, angle=90,hjust=0.5,align = TRUE, color="#39B600") +
  geom_cladelabel(node = 299, label = "italic(Cardiovirus)", parse=TRUE,offset=2.5, offset.text=0.4,  fontsize=4, angle=40,hjust=01,align = TRUE, color="#F8766D") +
  geom_cladelabel(node = 297, label = "italic(Mischivirus)",parse=TRUE,offset=2.5, offset.text=0.4,  fontsize=4,angle=14,hjust=0.8,align = TRUE, color="#00BF7F") +
  geom_cladelabel(node = 572, label = "italic(Teschovirus)",parse=TRUE,offset=2.5, offset.text=0.4,  fontsize=4,angle=356,hjust=0.3,align = TRUE, color="#E76BF3") +
  geom_cladelabel(node = 366, label = "italic(Hepatovirus)",parse=TRUE,offset=2.5, offset.text=0.4,  fontsize=4,angle=65,hjust=0.4,align = TRUE, color="#D89000") +
  geom_cladelabel(node = 328, label = "italic(Sapelovirus)",parse=TRUE,offset=2.5, offset.text=0.4,  fontsize=4,angle=330,hjust=0.3,align = TRUE, color="#00BFC4") +
  geom_cladelabel(node = 353, label = "italic(Bat_picornavirus)",parse=TRUE,offset=2.5, offset.text=0.4, angle=304,hjust=0.3, fontsize=4,align = TRUE, color="#9590FF") +
  geom_cladelabel(node = 463, label = "italic(Sapovirus)",parse=TRUE,offset=2.5, offset.text=0.4,angle=45, fontsize=4,align = TRUE,hjust=1, color="#00B0F6")
p1.1

##Add contig/read metadata
#attach various metadata to p1
p1.1 <- p1.1 %<+% region

#pop contig data on top of the tree
colz3 = c('Africa' =  "#0F6FC6", 'Asia' = "#009DD9", "Australia"="#0BD0D9","Europe"="#10CF9B",
          "North America"="#7CCA62","South America"="#A5C249")
p2<-p1.1+geom_fruit(#data=contig,
                  geom=geom_tile,
                  mapping=aes(fill=Region),
                  pwidth=1,
                  offset=0.15)+
                  scale_fill_manual(values=colz3)
                  #scale_fill_gradient(low="peachpuff1", high="orangered3")
                  
p2

#attach various metadata to p2
p2<-p2 %<+% host

colz4 = c('Bat' =  "#DE9ED6", 'Bird' = "#A55194", "Canine"="#7B4173","Elephant"="#CE6DBD",
          "Feline"="#E7969C","Hedgehog"="#D6616B","Human"="#AD494A","Lab strain"="black","Marsupial"="#843C39",
          "Non-human primate"="#E7CB94","Pig"="#E7BA52","Rodent"="#BD9E39","Shrew"="#7375B5","Ungulate"="#9C9EDE")
p3<-p2+new_scale_fill()+
  geom_fruit(#data=reads,
  geom=geom_tile,
  mapping=aes(fill=Host_class),
  pwidth=1,
  offset=0.15)+ 
  scale_fill_manual(values=colz4)+
  #scale_fill_gradient(low="lightblue1", high="royalblue2")+
  #scale_fill_viridis(option="B", name="Novel\nreads (log10)", direction = -1) +
 theme(
                 #legend.direction = "none",
                 legend.margin = margin(c(0,0,0,0)),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10),
                 legend.key.size = unit(0.3, "cm"),
       plot.margin = unit(c(0, 0, 0, 0), 
                          "cm"))
p3


#attach various metadata to p3
#p3<-p3 %<+% novel

base<-p3+new_scale_fill()+
  geom_fruit(#data=reads,
    geom=geom_tile,
    mapping=aes(fill=Seq_type),
    width=0.5,
    offset=0.1)  + 
  guides(fill = guide_legend(order=2))+
  scale_fill_manual(values=c("Novel seq"="tomato1","Reference seq"="grey88"), name="Seq type")+
  theme(
    #legend.direction = "none",
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = unit(c(0, 0, 0, 0), 
                       "cm"))
base


##Put the map and both summary figs together
library(cowplot)
library(ggplotify)

Fig1.1<-plot_grid(base, NULL, labels=c("",""),
                  rel_widths = c(1,1), rel_heights = c(3,0.5),
                  ncol=1, align="hv", axis="l", label_size = 23)
Fig1.1
Fig1.2<-as.ggplot(Fig1.2)

Fig1.2<-plot_grid(Fig1.1, sampling, labels=c("A","B"),
                  rel_widths = c(2.5,1.5), rel_heights = c(3,1),
                  ncol=2, align="hv", axis="l", label_size = 23)
Fig1.2
Fig1.2<-as.ggplot(Fig1.2)

ggsave(file = paste0(homewd, "/final_figures/Fig1_diversity_summary_phylogeny.pdf"),
       plot = Fig1.2,
       units="mm",  
       width=130, 
       height=70, 
       scale=3, 
       dpi=300)



