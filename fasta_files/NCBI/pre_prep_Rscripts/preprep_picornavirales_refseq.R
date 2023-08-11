rm(list=ls())

library(plyr)
library(dplyr)


#set wd
homewd = "/Users/gwenddolenkettenburg/Desktop/mada_bat_picornavirus"
setwd(paste0(homewd, "/fasta_files/NCBI/"))

#load the dataset and query
dat <- read.csv(file = "picornavirales_refseq_metadata.csv", header = T, stringsAsFactors = F)

head(dat)
#look at the unique hosts
sort(unique(dat$Host))
sort(unique(dat$Species[dat$Host==""])) #these ones have no hosts

#select all the bats
bat.dat = subset(dat, Host=="Chiroptera" | Host =="Hipposideros armiger" | Host =="Miniopterus magnater"|
                   Host=="Miniopterus pusillus" | Host == "Miniopterus schreibersii" | Host =="Myotis ricketti" | Host =="Pipistrellus pipistrellus"|
                   Host== "Rhinolophus sinicus" | Host== "Eidolon helvum" | Host== "Miniopterus" | Host=="Miniopterus fuliginosus")

sort(unique(bat.dat$Species))

nrow(bat.dat) 

#then select the ref seq for the hosts that are not bats
ref.dat <- subset(dat, Host=="") #check those with no host, have to manually fill in 

#now take just those with no bats...
ref.sub = subset(dat, Host!="Chiroptera" & Host!="Hipposideros armiger" & Host !="Miniopterus magnater" & Host !="Miniopterus pusillus" &
                   Host!="Miniopterus schreibersii" & Host!="Myotis ricketti" & Host !="Pipistrellus pipistrellus" & Host !="Rhinolophus sinicus" &
                   Host!="Eidolon helvum" & Host!="Miniopterus" & Host!="Miniopterus fuliginosus")

#now put these together to draw from GenBank
bat.picorna <- rbind(bat.dat, ref.sub)

#Check duplicated records
bat.picorna <- bat.picorna[!duplicated(bat.picorna),] #no duplicates

#and get the text to download from NCBI, check duplicated accession numbers
bat.picorna <- bat.picorna[!duplicated(bat.picorna$Accession),] #no duplicates


#and remove those that are repeats of the same record 
#all.bat.picorna = subset(all.CoV, Accession!="KU762338" &
                   # Accession!= "KF636752" &
                   # Accession!="GU190215" &
                   # Accession!="EF065505" &
                   # Accession!="EF065509"&
                   # Accession!="EF065513" &
                   # Accession!="KX574227")

accession_num <- paste(c(bat.picorna$Accession), collapse = ",")

#now put this into your webbrowser to download
text.for.NCBI <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=",accession_num)

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_015934.1,NC_015940.1,NC_015941.1,NC_017936.1,NC_025219.1,NC_028366.1,NC_030843.1,NC_031749.1,NC_033776.1,NC_033818.1,NC_033819.1,NC_033820.1,NC_033823.1,NC_033824.1,NC_034381.1,NC_038313.1,NC_038961.1,NC_043072.1,NC_000940.1,NC_001366.1,NC_001430.1,NC_001472.1,NC_001479.1,NC_001481.2,NC_001489.1,NC_001490.1,NC_001543.1,NC_001612.1,NC_001617.1,NC_001632.1,NC_001834.1,NC_001859.1,NC_001874.1,NC_001918.1,NC_001959.2,NC_002058.3,NC_002066.1,NC_002548.1,NC_002551.1,NC_002615.1,NC_003003.1,NC_003005.1,NC_003113.1,NC_003445.1,NC_003446.1,NC_003495.1,NC_003496.1,NC_003502.1,NC_003509.1,NC_003545.1,NC_003549.1,NC_003550.1,NC_003615.1,NC_003621.1,NC_003622.1,NC_003623.1,NC_003626.1,NC_003628.1,NC_003693.1,NC_003694.1,NC_003738.1,NC_003741.1,NC_003779.1,NC_003781.1,NC_003782.1,NC_003783.1,NC_003784.1,NC_003785.2,NC_003787.1,NC_003791.1,NC_003799.1,NC_003839.2,NC_003840.1,NC_003924.1,NC_003976.2,NC_003983.1,NC_003985.1,NC_003987.1,NC_003988.1,NC_003990.1,NC_004064.1,NC_004365.1,NC_004421.1,NC_004439.1,NC_004441.1,NC_004451.1,NC_004541.1,NC_004542.1,NC_004807.1,NC_004830.2,NC_005092.1,NC_005096.1,NC_005097.1,NC_005266.1,NC_005281.1,NC_005289.1,NC_005290.1,NC_006056.1,NC_006057.1,NC_006269.1,NC_006271.1,NC_006494.1,NC_006553.1,NC_006554.1,NC_006559.1,NC_006875.1,NC_006964.1,NC_006965.1,NC_007522.1,NC_007916.1,NC_008029.1,NC_008182.1,NC_008250.2,NC_008311.1,NC_008580.1,NC_008714.1,NC_009013.1,NC_009025.1,NC_009448.2,NC_009530.1,NC_009757.1,NC_009758.1,NC_009891.1,NC_009996.1,NC_010354.1,NC_010415.1,NC_010624.1,NC_010709.1,NC_010810.1,NC_010987.1,NC_011050.1,NC_011190.1,NC_011349.1,NC_011704.1,NC_012212.1,NC_012531.1,NC_012699.1,NC_012798.1,NC_012800.1,NC_012801.1,NC_012802.1,NC_012957.1,NC_012986.1,NC_013075.1,NC_013218.1,NC_013695.1,NC_014137.1,NC_014411.1,NC_014412.1,NC_014413.1,NC_014793.1,NC_015414.1,NC_015492.1,NC_015626.1,NC_015936.1,NC_016156.1,NC_016403.1,NC_016405.1,NC_016443.1,NC_016769.1
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_016964.1,NC_017939.1,NC_018226.1,NC_018383.1,NC_018384.1,NC_018400.1,NC_018506.1,NC_018570.2,NC_018613.1,NC_018668.1,NC_019712.1,NC_020898.1,NC_021178.1,NC_021201.1,NC_021220.1,NC_021482.1,NC_021566.1,NC_021567.1,NC_022004.1,NC_022332.1,NC_022611.1,NC_022798.1,NC_022802.1,NC_023016.1,NC_023021.1,NC_023022.1,NC_023162.1,NC_023422.1,NC_023483.1,NC_023627.1,NC_023676.1,NC_023858.1,NC_023861.1,NC_023984.1,NC_023985.1,NC_023987.1,NC_023988.1,NC_024016.1,NC_024031.1,NC_024070.1,NC_024073.1,NC_024078.1,NC_024120.1,NC_024489.1,NC_024497.1,NC_024765.1,NC_024766.1,NC_024767.1,NC_024768.1,NC_024769.1,NC_024770.1,NC_025114.1,NC_025432.1,NC_025474.1,NC_025479.2,NC_025675.1,NC_025676.1,NC_025788.1,NC_025835.1,NC_025890.1,NC_025961.1,NC_026249.1,NC_026250.1,NC_026315.1,NC_026316.1,NC_026470.1,NC_026733.1,NC_026921.1,NC_027026.1,NC_027054.1,NC_027122.1,NC_027128.1,NC_027214.1,NC_027713.1,NC_027818.1,NC_027915.1,NC_027917.1,NC_027918.1,NC_027919.1,NC_027926.1,NC_028139.1,NC_028363.1,NC_028364.1,NC_028365.1,NC_028380.1,NC_028964.1,NC_028970.1,NC_028981.1,NC_029038.1,NC_029052.1,NC_029306.1,NC_029307.1,NC_029309.1,NC_029645.1,NC_029646.1,NC_029647.1,NC_029854.1,NC_029905.1,NC_030454.1,NC_030651.1,NC_030793.1,NC_030886.1,NC_031105.1,NC_031106.1,NC_031324.1,NC_031338.1,NC_031687.1,NC_031688.1,NC_031766.1,NC_032087.1,NC_032112.1,NC_032115.1,NC_032126.1,NC_032222.1,NC_032270.1,NC_032978.1,NC_033081.1,NC_033152.1,NC_033492.1,NC_033619.1,NC_033695.1,NC_033793.1,NC_034206.1,NC_034217.1,NC_034245.1,NC_034267.1,NC_034384.1,NC_034385.1,NC_034444.1,NC_034453.1,NC_034617.1,NC_034971.1,NC_034973.1,NC_035110.1,NC_035115.1,NC_035184.1,NC_035198.1,NC_035214.1,NC_035218.1,NC_035221.1,NC_035450.1,NC_035455.1,NC_035457.1,NC_035675.1,NC_035779.1,NC_036389.1,NC_036585.1,NC_036588.1,NC_037654.1,NC_038301.1,NC_038302.1,NC_038303.1,NC_038304.1,NC_038305.1,NC_038306.1,NC_038307.1,NC_038308.1,NC_038309.1,NC_038310.1,NC_038311.1,NC_038312.1,NC_038314.1,NC_038315.1,NC_038316.1,NC_038317.1,NC_038318.1,NC_038319.1,NC_038320.1,NC_038321.1,NC_038760.1,NC_038761.1,NC_038762.1,NC_038764.1,NC_038765.1,NC_038767.1,NC_038878.1,NC_038880.1,NC_038957.1
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_038989.1,NC_039004.1,NC_039073.1,NC_039077.1,NC_039209.1,NC_039210.1,NC_039212.1,NC_039235.1,NC_039236.1,NC_039475.1,NC_039476.1,NC_039477.1,NC_039897.1,NC_040399.1,NC_040417.1,NC_040574.1,NC_040586.1,NC_040587.1,NC_040601.1,NC_040605.1,NC_040611.1,NC_040642.1,NC_040646.1,NC_040673.1,NC_040674.1,NC_040675.1,NC_040684.1,NC_040716.1,NC_040724.1,NC_040832.1,NC_040876.1,NC_043447.1,NC_043512.1,NC_043516.1,NC_043519.1,NC_043542.1,NC_043544.1,NC_043684.1,NC_044045.1,NC_044046.1,NC_044047.1,NC_044853.1,NC_044854.1,NC_044855.1,NC_044856.1,NC_044932.1,NC_045762.1,NC_055108.1,NC_055125.1,NC_055156.1,NC_055159.1,NC_055160.1,NC_055161.1,OP287812.1,NC_040798.1,NC_011829.1


#once downloaded, send to MAFFT for alignment
rm(list=ls())
#then, after alignment is ready, prepare the names for RAxML (no space, semicolon, colon, parentheses, dash, slash, comma, quote allowed in name (should just all be underscore)
library(seqinr)
#library(msa)
alignment1 <- read.alignment(file = "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/fasta_files/NCBI/trim_align_picornavirales.fasta", format="fasta", forceToLower = F)


tmp <- as.list(alignment1$nam)

change.spacing <- function(df){
  df_new <- sapply(strsplit(df,"-"), function(x) x[[1]])
  return(df_new)
}

names_new = c(unlist(lapply(tmp, change.spacing)))

#new_names <- sub("__", "_", new_names) 
class(alignment1$seq)
write.fasta(sequences = as.list(alignment1$seq), names = names_new, file.out =  "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/fasta_files/NCBI/cleaned_picornavirales.fasta", as.string = T, open="w")

#now send to modeltest and eventually RAxML
