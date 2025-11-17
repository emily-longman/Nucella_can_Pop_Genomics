#### Moments analysis for SWD

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load data
SS_3pop <- list.files(path="data/processed/pop_structure/o3step/", pattern =  "*_output.3step.txt")
SS_3pop_v = as.vector(unlist(lapply(SS_3pop, function(x) paste0('data/processed/pop_structure/o3step/', x))))


Admix_3pop <- list.files(path="data/processed/pop_structure/admix/", pattern =  "*_output.admix.txt")
Admix_3pop_v = as.vector(unlist(lapply(Admix_3pop, function(x) paste0('data/processed/pop_structure/admix/', x))))


pop3models_SS = foreach(i=SS_3pop, .combine = "rbind")%do%{
  tmp <- fread(i)
  
  data.table(Pair_name = tmp$Pair_name,
             AIC = tmp$AIC) %>%
    mutate(model = "SS")
}

pop3models_Admix = foreach(i=Admix_3pop, .combine = "rbind")%do%{
  tmp <- fread(i)
  
  data.table(Pair_name = tmp$Pair_name,
             AIC = tmp$AIC) %>%
    mutate(model = "Admix")
}

# ================================================================================== #

# Combine data
rbind(pop3models_SS, pop3models_Admix) ->
  AIC_3pop_models

# ================================================================================== #

AIC_3pop_models %>%
  group_by(model,Pair_name) %>%
  summarise(m = mean(log10(AIC)),
            sd = sd(log10(AIC)))

AIC_3pop_models %>%
  group_by(model,Pair_name) %>%
  slice_min(AIC)


AIC_3pop_models %>%
  group_by(model,Pair_name) %>%
  summarise(m = mean(log10(AIC)),
            sd = sd(log10(AIC))) %>%
  ggplot(aes(
    x=Pair_name,
    y=m,
    ymin = m - sd,
    ymax = m + sd,
    color = model
  )) + geom_errorbar(width = 0.08,  position=position_dodge(width=0.5)) + 
  geom_point(size = 2,  position=position_dodge(width=0.5)) + theme_bw() ->
  plot3pop_models

ggsave(plot3pop_models, 
       w=4, h =4,
       file = "output/figures/plot3pop_models.pdf")

AIC_3pop_models %>%
  group_by(model,Pair_name) %>%
  slice_min(AIC)

# ================================================================================== #  
# ================================================================================== #

# Extract estimates for admixed model

# Load all data
all_Admix = foreach(i=Admix_3pop_v, .combine = "rbind")%do%{
  tmp <- fread(i)
  
  tmp %>%
    mutate(model = "Admix")
}

# ================================================================================== #

mu=2.8e-9
g=1

all_Admix %>%
  group_by(Pair_name) %>%
  slice_min(AIC) %>%
   as.data.frame() %>%
  mutate(Nref= theta/(4*mu*1585293896)) %>%
  mutate(divergence_time = 2*Nref*T_split*g,
         admix_time = 2*Nref*T_admix*g,
         pop1size = nu1*(theta/(4*mu*L)),
         pop2size = nu2*(theta/(4*mu*L)),
         pop3size = nu3*(theta/(4*mu*L)),
         ) 

# Write table
write.csv(all_Admix, "output/tables/moments_all_Admix.csv")