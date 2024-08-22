library(Seurat)
library(drake)
library(tidyverse)
library(cowplot)
library(data.table)
library(tradeSeq)
library(ROGUE)
library(pals)
library(readxl)
library(harmony)
library(tidyverse)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Nebulosa)
library(scProportionTest)
library(presto)
library(fgsea)
library(msigdbr)
library(ReactomePA)
library(biomaRt)
library(BioNERO)

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)




newpan <- c("#DBACAC","#9593BC","#80B379","#2C8020","#F7C070")

StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          colors = tableau_color_pal(palette = "Tableau 20")(20),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x,pt.size = pt.size, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.ticks.x = element_line(),axis.text.x = element_text(angle = 90))
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y)+scale_fill_manual(values = colors))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
