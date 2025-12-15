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
SS_3pop <- list.files(path="data/processed/pop_structure/o3step_relaxed_filt/", pattern =  "*_output.3step.txt")
SS_3pop_v = as.vector(unlist(lapply(SS_3pop, function(x) paste0('data/processed/pop_structure/o3step_relaxed_filt/', x))))


Admix_3pop <- list.files(path="data/processed/pop_structure/admix_relaxed_filt/", pattern =  "*_output.admix.txt")
Admix_3pop_v = as.vector(unlist(lapply(Admix_3pop, function(x) paste0('data/processed/pop_structure/admix_relaxed_filt/', x))))

pop3models_SS = foreach(i=SS_3pop_v, .combine = "rbind")%do%{
  tmp <- fread(i)
  
  data.table(Pair_name = tmp$Pair_name,
             AIC = tmp$AIC) %>%
    mutate(model = "SS")
}

pop3models_Admix = foreach(i=Admix_3pop_v, .combine = "rbind")%do%{
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

pop3models_Admix %>%
  group_by(model,Pair_name) %>%
  summarise(m = mean(log10(AIC)),
            sd = sd(log10(AIC)))

AIC_3pop_models %>%
  group_by(model,Pair_name) %>%
  summarise(m = mean(log10(AIC)),
            sd = sd(log10(AIC)))

AIC_3pop_models %>%
  group_by(model,Pair_name) %>%
  summarise(m = mean(AIC),
            sd = sd(AIC))

pop3models_Admix %>%
  group_by(model,Pair_name) %>%
  slice_min(AIC)


colors = c("darkorange2",  "mediumorchid4")

# Graph 
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
  geom_point(size = 2,  position=position_dodge(width=0.5)) + 
  scale_color_manual(values=colors) + 
  xlab(expression(paste(italic("N. canaliculata"), " Population"))) + ylab(expression(Log[10](AIC))) + theme_bw() ->
  plot3pop_models

pdf("output/figures/pop_structure/plot3pop_models_relaxed_filt.pdf", width = 6, height = 5)
plot3pop_models
dev.off()

AIC_3pop_models %>%
  group_by(model,Pair_name) %>%
  slice_min(AIC)
  
# ================================================================================== #
# ================================================================================== #

# Extract estimates for admixed model
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
write.csv(all_Admix, "output/tables/moments_all_Admix_relaxed_filt.csv")