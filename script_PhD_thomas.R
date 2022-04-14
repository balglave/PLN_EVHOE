library(FactoMineR)

# Shape data
load("data/PhD_Outrequin/invertebrate.Rdata")

sampling <- read.csv(file = "data/PhD_Outrequin/evhoeALL_TRAWLRESUM.txt",header=T,sep="\t")

yrs<-seq(2008,2020)
spe<-invertebrate %>%
  filter(position=="epi") %>%
  filter(year%in%yrs) %>%  
  dplyr::select(unique_id,scientific_name,abundance,sampling_point) %>%
  # left_join(sampling[,c('unique_id','surface_corr')])%>%
  filter(!is.na(abundance)) # %>%
  # filter(abundance>0)

sum_spe <- spe %>%
  group_by(scientific_name) %>%
  dplyr::summarise(abundance = sum(abundance))

spe_2 <- spe %>%
  filter(!scientific_name %in% sum_spe$scientific_name[which(sum_spe$abundance==0)]) %>%
  pivot_wider(names_from = scientific_name,
              values_from = abundance)

spe_2[is.na(spe_2)] <- 0
spe_2 <- spe_2[-329,] # no positive value
Abundance_2 <- spe_2[,-c(1,2)]

# Plot correlation
cor_df <- cor(log(Abundance_2[,]+1))
corrplot::corrplot(corr = cor_df,
                   addgrid.col = NA,
                   tl.cex=0.65,
                   cl.pos="n",
                   cl.lim=c(-1,1))

model_inputs <- prepare_data(counts = Abundance_2,
                             covariates = rep(1,nrow(Abundance_2)))

# Fit model
res <- PLNPCA(Abundance ~ 1,data=model_inputs,ranks=1:7)
plot(res)
res_best <-  getBestModel(res,"BIC")
plot(res_best,nb_axes = 2)

# save(file="res/PhD_Outrequin/res_best.RData",data=res_best)
load("res/PhD_Outrequin/res_best.RData")

# Make HCPC
PC_df <- data.frame(PC=res_best$scores,unique_id=spe_2$unique_id)
# save(file="res/PhD_Outrequin/PC_df.RData",data=PC_df)
res <- HCPC(as.data.frame(res_best$scores))

