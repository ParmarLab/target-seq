figure3 <- function(outdir,hpsc.thneuron.file="processed_data/hpsc.human_harmony.THneurons.rds"){
  
  dir.create(outdir)
  dir.create(paste0(outdir,"/figure5/"))
  outdir <- paste0(outdir,"/figure5/")
  
  # Prepare hpsc data
  hpsc.thneuron.data <- readRDS(hpsc.thneuron.file)
  hpsc.thneuron.data$NewCluster <- factor(paste0("DA-",as.integer(hpsc.thneuron.data$seurat_clusters)))
  hpsc.thneuron.data$cellbc_sample <- colnames(hpsc.thneuron.data)
  
  
  # Figure 5A
  # Export data for spring
  system("mkdir -p spring")
  data_to_write_out <- as.matrix(GetAssayData(hpsc.thneuron.data[,],assay =  "RNA"))
  # expression
  fwrite(x = data_to_write_out, file = "spring/expression.csv.gz",compress = "gzip")
  #genelist
  fwrite(x = data.frame(rownames(hpsc.thneuron.data)), 
         file = "spring/genelist.csv.gz",compress = "gzip",row.names = F,col.names = F)
  # meta
  groups<-t(data.frame(c("Cluster",new_df$NewCluster),
                       c("Age",as.character(new_df$Age))))
  fwrite(x =groups , file = "spring/day.groups.csv",row.names = F,col.names = F,quote = F)
  
  
  # Figure 5B import from spring and run slingshot
  
  coord <- read.csv("spring/coordinates.txt",header = F)
  edge <- read.csv("spring/edge_list.txt",header = F)
  sd <- read.csv("spring/day.groups.csv")
  sd <- t(sd)[-1,]
  colnames(sd)<-c("Location" ,"NamedClusters" ,"Sample")
  coord <- cbind(coord,sd)
  rd <- coord[,c("V2","V3")]
  rd$V2 <- rd$V2*-1
  cl <- coord$NamedClusters
  pto <- slingshot(rd, cl, start.clus = '0')
  sds <- as.SlingshotDataSet(pto)
  plot(rd, col = tableau_color_pal()(10)[factor(cl)], asp = 1,pch=19)
  lines(sds, type = 'l', lwd = 3,show.constraints=T)

  # Figure 5C
  counts <- read.csv("spring/expression.csv.gz")
  genes <- read.csv("spring/genelist.csv.gz",header = F)
  rownames(counts)<-genes$V1
  filt_counts <- as.matrix(counts[rowSums(counts > 1), ])
  sce <- fitGAM(counts = as.matrix(filt_counts), sds = sds,)
  
  startRes <- startVsEndTest(sce,lineages=T)
  
  topGenes <- c("TMTC2","PLXDC2","MGAT4C","KCNIP4","LUZP2","SGCZ","KCNH8","DPP10")
  
  plot_grid(plotGeneCount2(sds, filt_counts,  gene=topGenes[1], models = sce,which_curves = 3)+ggtitle(topGenes[1])+scale_colour_viridis()&NoLegend(),
            plotGeneCount2(sds, filt_counts,  gene=topGenes[2], models = sce,which_curves = 3)+ggtitle(topGenes[2])+scale_colour_viridis()&NoLegend(),
            plotGeneCount2(sds, filt_counts,  gene=topGenes[3], models = sce,which_curves = 3)+ggtitle(topGenes[3])+scale_colour_viridis()&NoLegend(),
            plotGeneCount2(sds, filt_counts,  gene=topGenes[4], models = sce,which_curves = 3)+ggtitle(topGenes[4])+scale_colour_viridis()&NoLegend(),
            plotGeneCount2(sds, filt_counts,  gene=topGenes[5], models = sce,which_curves = 3)+ggtitle(topGenes[5])+scale_colour_viridis()&NoLegend(),
            plotGeneCount2(sds, filt_counts,  gene=topGenes[6], models = sce,which_curves = 3)+ggtitle(topGenes[6])+scale_colour_viridis()&NoLegend(),
            plotGeneCount2(sds, filt_counts,  gene=topGenes[7], models = sce,which_curves = 3)+ggtitle(topGenes[7])+scale_colour_viridis()&NoLegend(),
            plotGeneCount2(sds, filt_counts,  gene=topGenes[8], models = sce,which_curves = 3)+ggtitle(topGenes[8])+scale_colour_viridis()&NoLegend(),
            ncol = 4)
  
  # Figure 5D
  
  grn.a9 <- read.csv("processed_data/raw_GRN_A9.csv")
  grn.a9.fil <- grn.a9 %>% filter(X.logp > 14)
  
  # most central
  df1 <- data.frame(table(grn.a9.fil$source))
  df1 <- df1[df1$Freq>25,]
  grn.a9.fil2 <- grn.a9.fil[grn.a9.fil$source %in% df1$Var1,]
  
  plot_grn(grn.a9.fil2[,-c(1)],show_labels = "tophubs",top_n_hubs = 35,ranked = F)+
    scale_color_manual(values = c("darkgrey","#D3335C"))+ggtitle("A9")


  # Figure 5E
  
  grn.a10 <- read.csv("processed_data/raw_GRN_A10.csv")
  grn.a10.fil <- grn.a10 %>% filter(X.logp > 13)
  
  # most central
  df1 <- data.frame(table(grn.a10.fil$source))
  df1 <- df1[df1$Freq>25,]
  grn.a10.fil2 <- grn.a10.fil[grn.a10.fil$source %in% df1$Var1,]
  
  plot_grn(grn.a10.fil2[,-c(1)],show_labels = "tophubs",top_n_hubs = 35,ranked = F)+
    scale_color_manual(values = c("darkgrey","#669DD1"))+ggtitle("a10")
  
  
 
  # Figure 5F
  
  # write for celloracle
  coord <- coord[,-1]
  colnames(coord)[c(1,2)]<-c("UMAP1","UMAP2")
  rownames(coord)<-coord$Cell
  data.table::fwrite(coord,file = "metadata.co.csv")
  data.table::fwrite(counts,file = "counts.co.csv")
  # celloracle ran in separate notebook
  
  # Figure 5I
  
  
  
  
  
  }
