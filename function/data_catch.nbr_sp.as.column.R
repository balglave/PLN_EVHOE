#####################################################################
## Modify survey counts data so that count per species are in columns
#####################################################################

data_catch.nbr_sp.as.column <- function(data_survey.nbr){
  
  
  ## delete duplicated data
  data_survey.nbr <- data_survey.nbr[!duplicated(data_survey.nbr),]
  data_survey.nbr %>% 
    dplyr::select("Survey",
                  "StNo",
                  "HaulNo",
                  "Year",
                  "Stratum",
                  "Distance",
                  "HaulNo",
                  "scientificname",
                  "Sex",
                  "CatIdentifier",
                  "TotalNo",
                  "lati",
                  "long") %>%
    mutate(scientificname = paste("TotalNo_",scientificname,sep='')) -> tmp
  
  # create toto to save caracteristics of the hauls
  tmp %>% dplyr::select(-c(TotalNo,scientificname,CatIdentifier,Sex)) -> toto
  
  # take only 1st data row for each haul : others are duplicates
  toto %>% group_by(HaulNo,Year) %>% slice(1) -> toto
  
  # species as columns
  tmp %>%
    group_by(scientificname,HaulNo,Year) %>%
    dplyr :: summarise(TotalNo = sum(as.numeric(TotalNo))) %>%
    spread(key= scientificname, value = TotalNo) -> tmp
  tmp <- full_join(toto,tmp, by = c("HaulNo","Year"))
  tmp[is.na(tmp)] <- 0
  return(tmp)
}
