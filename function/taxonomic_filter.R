taxinomic_filter<-function(catch_data,pelagic_rm=TRUE,cephalopoda=FALSE, poissons=FALSE,invertebrate=FALSE){
  catch<-catch_data
  cat<-catch
  file_taxo_path<-"C:/Users/toutrequ/AOBox/thèse/AFH/project/function/taxonomic_filter/"
  load(paste0(file_taxo_path,"DATAtaxo.Rdata"))
  load(paste0(file_taxo_path,"DATAtaxoMISS.Rdata"))
  load(paste0(file_taxo_path,"DATAtaxoDIFF.Rdata"))
  taxo <- DATAtaxo %>% select(AphiaID, Nom_Scientifique,genus,family,order,class,phylum) %>% distinct()
  cat<-cat %>% left_join(taxo,by=c("scientific_name"="Nom_Scientifique"))
  
  load(paste0(file_taxo_path,"to_remove.Rdata"))
  
  # cat <-cat %>% filter(!order%in%c("Asparagales","Caryophyllales"))
  # pelagos<-c("Semaeostomeae","Pteropoda","Rhizostomeae","Coronatae",
  #            "Pyrosomatida","Narcomedusae","Cydippida","Argentiniformes" ,
  #            "Scombriformes","Carangiformes","Salpida","Euphausiacea",
  #            "Mysidacea","Clupeiformes","Atheriniformes","Lamniformes",
  #            "Beloniformes","Acanthuriformes")
  cat<-cat %>% filter(!scientific_name%in%to_rm$scientific_name)
  
  
  if (pelagic_rm){
    load(paste0(file_taxo_path,"pelagique.Rdata"))
    cat<-cat %>% filter(!scientific_name%in%pelagique$scientific_name) 
  }
  ### tri pélagique et incohérences
  
  if(invertebrate){
    load(paste0(file_taxo_path,"invertébrés.Rdata"))
    cat<-cat %>% filter(scientific_name%in%invertébrés$scientific_name)
    return(cat)
  }
  
  if(poissons){
    load(paste0(file_taxo_path,"poissons.Rdata"))
    cat<-cat %>% filter(scientific_name%in%poissons$scientific_name)

    return(cat)
  }
  if(cephalopoda){
    load(paste0(file_taxo_path,"cephalo.Rdata"))
    cat<-cat %>% filter(scientific_name%in%cephalo$scientific_name)
    return(cat)
  }
  return(cat)
}



# drv = dbDriver("PostgreSQL")
# ben <- dbConnect(drv, host="sirs.agrocampus-ouest.fr", user="outrequin", password="MqaMsGx", dbname="bentroph")
# catch<-dbGetQuery(ben, "select * from backup.catch")
# esp_tot<-catch %>% filter(survey=="EVHOE") %>% dplyr::select(scientific_name,valid_name,aphia_id) %>% distinct()
# dim(esp_tot)
# esp_tot<-esp_tot %>% left_join(taxo,by=c("scientific_name"="Nom_Scientifique"))
# setwd("C:/Users/toutrequ/AOBox/thèse/AFH/project/data/taxo/EVHOE")
# # a retirer
# 
# to_rm<-esp_tot %>% filter(order%in%c("Asparagales","Caryophyllales")|scientific_name%in%c("Plantae","Hirudinea",
#                                                                                           "Animalia"))
# 
# # save(to_rm, file = "./to_remove.Rdata")
# # write.csv2(to_rm, file = "./format_csv/to_remove.csv",row.names = F)
# ### pelagique
# 
# 
# pelagos<-c("Semaeostomeae","Pteropoda","Rhizostomeae","Coronatae",
#            "Pyrosomatida","Narcomedusae","Cydippida","Argentiniformes" ,
#            "Scombriformes","Carangiformes","Salpida","Euphausiacea",
#            "Mysidacea","Clupeiformes","Atheriniformes","Lamniformes",
#            "Beloniformes","Acanthuriformes")
# pelagique<-esp_tot %>% filter(order%in%c('Myopsida','Oegopsida')|str_detect(scientific_name, "glauca")|str_detect(scientific_name, "Galeorhinus")|str_detect(scientific_name, "sedentaria")|order%in%pelagos|scientific_name%in%c("Scyphozoa")|str_detect(scientific_name,"equorea"))
# # save(pelagique, file = "./pelagique.Rdata")
# # write.csv2(pelagique, file = "./format_csv/pelagique.csv",row.names = F)
# #### poissons
# cati<-esp_tot %>% filter(order%in%c("Pleuronectiformes","Petromyzontiformes","Myxiniformes","Callionymiformes","Acanthuriformes",
#                                 "Perciformes","Notacanthiformes" , "Salmoniformes","Ateleopodiformes",
#                                 "Anguilliformes" , "Rajiformes","Mugiliformes","Tetraodontiformes",
#                                 "Torpediniformes","Syngnathiformes",
#                                 "Ophidiiformes","Squaliformes","Beryciformes","Chimaeriformes","Myliobatiformes","Blenniiformes",
#                                 "Stomiiformes","Alepocephaliformes","Hexanchiformes","Trachichthyiformes",
#                                 "Acropomatiformes","Aulopiformes",  "Lophiiformes" ,"Gadiformes","Carcharhiniformes",
#                                 "Zeiformes","Eupercaria incertae sedis",
#                                 "Mulliformes","Gobiiformes","Myctophiformes"),
#                                            !scientific_name%in%c("Animalia","Cephalopoda","Pisces","Actinopterygii"),!is.na(scientific_name)) %>% 
#   filter(!scientific_name%in%to_rm$scientific_name,!order%in%pelagique$scientific_name,!scientific_name%in%pelagique$scientific_name)
# cata<-esp_tot %>% filter(scientific_name%in%c("Pisces","Actinopterygii"))
# poissons<-rbind(cati,cata)
# # save(poissons, file = "./poissons.Rdata")
# # write.csv2(poissons, file = "./format_csv/poisson.csv",row.names = F)
# # ### cephalo
# cati<-esp_tot %>% filter(order%in%c("Sepiida","Octopoda"),
#                      !scientific_name%in%c("Scyphozoa","Plantae","Hirudinea",
#                                            "Animalia","Cephalopoda"),!is.na(scientific_name))
# cata<-esp_tot %>% filter(scientific_name%in%c("Cephalopoda"))
# cephalo<-rbind(cati,cata)
# # save(cephalo, file = "./cephalo.Rdata")
# # write.csv2(cephalo, file = "./format_csv/cephalo.csv",row.names = F)
# # ####invert
# 
# cat<-esp_tot %>% filter(!order%in%c("Pleuronectiformes","Petromyzontiformes","Myxiniformes","Callionymiformes","Acanthuriformes",
#                                 "Perciformes","Octopoda","Argentiniformes","Notacanthiformes" , "Salmoniformes","Ateleopodiformes",
#                                 "Anguilliformes" , "Clupeiformes","Rajiformes","Mugiliformes","Tetraodontiformes","Atheriniformes",
#                                 "Torpediniformes","Syngnathiformes",
#                                 "Ophidiiformes","Squaliformes","Beryciformes","Chimaeriformes","Myliobatiformes","Blenniiformes",
#                                 "Lamniformes","Beloniformes","Stomiiformes","Alepocephaliformes","Hexanchiformes","Trachichthyiformes",
#                                 "Acropomatiformes","Aulopiformes",  "Lophiiformes" ,"Gadiformes","Carcharhiniformes",
#                                 "Sepiida","Zeiformes","Eupercaria incertae sedis","Scombriformes","Carangiformes",
#                                 "Mulliformes","Myopsida" ,"Gobiiformes","Myctophiformes","Oegopsida",NA)&
#                     !scientific_name%in%c("Actinopterygii","Scyphozoa","Plantae","Hirudinea",
#                                           "Animalia","Cephalopoda","Pisces")&!is.na(scientific_name)&
#                     !str_detect(scientific_name,"equorea"))
# cat<-cat %>% filter(!scientific_name%in%to_rm$scientific_name&!scientific_name%in%cephalo$scientific_name&!scientific_name%in%pelagique$scientific_name)
# invertébrés<-cat
# # save(invertébrés, file = "./invertébrés.Rdata")
# # write.csv2(invertébrés, file = "./format_csv/invertébrés.csv",row.names = F)
