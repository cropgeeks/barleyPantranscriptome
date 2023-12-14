library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(readxl)
library(ggthemes)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(stringr)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

# read data
dat <- read_excel(args[1], sheet = 'Sheet1')
data2 <- read_excel(args[2], sheet = 'Sheet1')
#data2$gene_Set=gsub("ALL_GENES (1,5001 genes) corr= 0.395","ALL_GENES (15,001 genes) corr= 0.395",x =data2$gene_Set,fixed = T )
dat2=read_excel(args[3],'Sheet1')
dat2=dat2[c(1,2,3,11,12,13)]
dat_melt=melt(dat2, id.vars=c("COEF_VAR__FRACTION"))
tr=as.character(sort(unique(dat_melt$variable)))
dat_melt$variable=factor(dat_melt$variable,levels = tr)
dat_melt$line=ifelse(dat_melt$variable==tr[1],2,1)

# Filter data
geno <- c('Akashinriki', 'B1K.04.12', 'GoldenPromise', 'HOR3081', 'RGTPlanet')
dat2 <- dat %>% filter(GENOTYPE %in% geno, CONDITION == "ALL_GENES", REGION != "F")
#define font size
lab <- 5

# plot1
dat2$REGION=paste("REGION: ",dat2$REGION)
p1 <- ggplot(dat2, aes(x = as.numeric(as.factor(CIS_ELEM_name)), y = percent_INTERSECT_perf_hits,color= GENOTYPE)) +
  geom_line() +
  facet_wrap(~REGION, ncol = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        legend.position = c(0.9,0.3),
        strip.text = element_text(size = lab),
        axis.text = element_text(size = lab,margin=margin(b=-10)),
        axis.title = element_text(size = lab),
        legend.text = element_text(size = lab),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.background = element_rect(fill = 'white',color='white'),
        panel.background = element_rect(fill = 'whitesmoke'),
        panel.grid.major = element_blank(),
        panel.grid.minor =element_blank()) +
  scale_y_continuous(limits = c(10, 100))+
  scale_x_continuous(breaks=seq(1,30,2))+
  labs(x = NULL, y = 'Percentage (%)') +
  guides(color = guide_legend(title = NULL,byrow = TRUE)) +
  scale_color_manual(values = c('Akashinriki' = 'red', 'B1K.04.12' = 'blue', 'GoldenPromise' = 'green', 'HOR3081' = 'purple', 'RGTPlanet' = 'orange'))

# plot2
p2 <- ggplot(data2, aes(x = percent_upstream_TFBS_identity/100, y = percent_coherent_expression_pairs/100, color = str_wrap(gene_Set,20))) +
  geom_point(size = 1) +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_minimal() +
  theme(axis.text = element_text(size = lab),
        axis.title = element_text(size = lab),
        legend.text = element_text(size = lab),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor =element_blank(),
        panel.background = element_rect(fill = 'whitesmoke')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.25, 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.3, 0.65)) +
  labs(x = 'Percentage of TFBS identity (%)', y = 'Percentage of coherent expression (%)') +
  guides(color = guide_legend(title = NULL))

#plot 3
p3 <- ggplot(dat_melt, aes(x = as.numeric(factor(COEF_VAR__FRACTION)), y = value/100,
                           color = str_wrap(variable,width = 20))) +
  geom_line(aes(size=str_wrap(variable,width = 20)),show.legend = T,inherit.aes = T)+guides(size='none')+
  geom_point(aes(x = as.numeric(factor(COEF_VAR__FRACTION)), y = value/100,
                  color = str_wrap(variable,width = 20))) + 
  scale_size_manual(values = c(1,1,1,3,1))+
  scale_fill_discrete(breaks=tr)+
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = lab),
        legend.text = element_text(size = lab),
        legend.title = element_blank(),
        #legend.position = c(1,0.3),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'gray', size = 0.2)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(1:5),labels = c("Gene set 1","Gene set 2","Gene set 3","Gene set 4","Gene set 5"))+
  labs(x = 'Gene sets showing increasing coefficient of variation',
       y = 'Percentage (%)') +
 
  guides(color = guide_legend(title = NULL)) +
  scale_linetype_manual(values = c("#3C5488B2", "#4DBBD5B2", "#00A087B2", "#E64B35B2", "#F39B7FB2", "#8491B4B2","#91D1C2B2", "#DC0000B2", "#7E6148B2"))+
  scale_color_manual(values = c("#3C5488B2", "#4DBBD5B2", "#00A087B2", "#E64B35B2", "#F39B7FB2", "#8491B4B2","#91D1C2B2", "#DC0000B2", "#7E6148B2"))

#combining plots
z2=grid.arrange(p1,                     
             p2, p3,                               
             ncol = 5, nrow = 5, 
             layout_matrix = rbind(c(1,1,1,1,1),c(1,1,1,1,1),c(1,1,1,1,1),c(2,2,3,3,3), c(2,2,3,3,3)))
all_c <- as_ggplot(z2) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 5,
                  x = c(0, 0, 0.42), y = c(1, 0.42, 0.42))
cm = 1/2.54
# save image with combined plots
ggsave(filename = args[4], plot =all_c, width = 18 * cm, height = 11.5 * cm, dpi = 900, 
       bg = "white",device = "png",units = "in",limitsize = F)
