library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(poolfstat)
library(FactoMineR)
require(foreach)
require(gmodels)

files <- system(paste("ls slim_glacial_refugiaNorth" ), intern = T)

inf <- foreach(i=files, .combine = "rbind", .errorhandling = "remove")%do%{
  tmp <- fread(paste("slim_glacial_refugiaNorth/",i, sep = ""))
}
inf %<>% 
  mutate(model = "NorthRef")
names(inf) = c(paste("p", 0:6, sep = "" ), "admix", "model")


files2 <- system(paste("ls slim_glacial_refugiaSouthExpansion" ), intern = T)

inf2 <- foreach(i=files2, .combine = "rbind", .errorhandling = "remove")%do%{
  tmp <- fread(paste("slim_glacial_refugiaSouthExpansion/",i, sep = ""))
}
inf2 %<>% 
  mutate(model = "SouthExp")
names(inf2) = c(paste("p", 0:6, sep = "" ), "admix", "model")


rbind(inf, inf2) %>%
  melt(id = c("admix", "model")) %>%
  separate(variable, remove = F, 
           into = c("p","simlat"), sep = 1) %>%
  ggplot(aes(
    x=variable,
    y=value,
    fill = model
  )) + geom_boxplot() + theme_bw() ->
  het_plot_box

ggsave(het_plot_box, file = "het_plot_box.pdf",
       h=4, w=5)

