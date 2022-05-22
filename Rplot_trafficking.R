library(readxl)
library(ggplot2)
library(magrittr)
library(ggpubr)
results <- read_excel("X:/91 Data and analysis/YJunqi/Sensitivity/20190103 HEK293T cells TMR(Halo tag) for Voltron_D81S +2APB/5uM 15min labeling Dish2/cell4/170448_CQ/Trafficking.xlsx")
results['Y'] = results['Y']/max(results['Y'])

traffic_plot = ggline(results,'X','Y',plot_type = c("l"),size = 1)+ylab('Normanized intensity')+xlab('Distance(pixels)')+ scale_x_discrete(breaks = seq(2,nrow(results), 10))+
  theme(axis.title.x = element_text(size = 16,family = "arial"))+
  theme(axis.text.x = element_text(size = 16, family = "arial"))+
  theme(axis.title.y = element_text(size = 16, family = "arial"))+
  theme(axis.text.y = element_text(size = 16, family = "arial")