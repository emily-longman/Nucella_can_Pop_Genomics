#### Moments analysis for SWD

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)


SS_3pop <- system("ls o3step/*_output.3step.txt", intern = T)
Admix_3pop <- system("ls admix/*_output.admix.txt", intern = T)

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

rbind(pop3models_SS, pop3models_Admix) ->
  AIC_3pop_models

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
       file = "plot3pop_models.pdf")

AIC_3pop_models %>%
  group_by(model,Pair_name) %>%
  slice_min(AIC)
  
####
all_Admix = foreach(i=Admix_3pop, .combine = "rbind")%do%{
  tmp <- fread(i)
  
  tmp %>%
    mutate(model = "Admix")
}
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
