figure2 <- function(outdir,hpsc.neuron.file,fetal.neuron.file) {
  
  dir.create(outdir)
  dir.create(paste0(outdir,"/figure2/"))
  outdir <- paste0(outdir,"/figure2/")
  
  # Prepare hpsc data
  hpsc.neuron.data <- readRDS(hpsc.neuron.file)
  hpsc.neuron.data$NewCluster <- as.integer(hpsc.neuron.data$seurat_clusters)
  Idents(hpsc.neuron.data) <- hpsc.neuron.data$NewCluster
  
  hpsc.thneuron.data <- readRDS("processed_data/hpsc.human_harmony.THneurons.rds")
  hpsc.thneuron.data$NewCluster <- factor(paste0("DA-",as.integer(hpsc.thneuron.data$seurat_clusters)))
    
  # Figure 2A
  DimPlot(hpsc.neuron.data,label = T,raster = F,pt.size = 1)+scale_color_manual(values = stevens.greenblue())
  ggsave(paste0(outdir,"1A.umap.png"),w=8,h=6)
  
  
  # Figure 2B
  genelist2 <- data.frame(read_excel("genelists//genelist_violin_fig2.xlsx"))
  
  for(c in unique(genelist2$Category)){
    
    StackedVlnPlot(hpsc.neuron.data,
                   features =genelist2[genelist2$Category==c,"Gene"],
                   colors = stevens.greenblue())&scale_fill_manual(values = stevens.greenblue())
    ggsave(paste0(outdir,"/1B.vln.neurons.",c,".pdf"),w=8,h=8)
  }
  
  # Figure 2C
  
  fetal<-readRDS(fetal.file)
  fetal <- RenameIdents(fetal, "0"="Astro","1"="FPP", "2"="Neuron","3"="Astro","4"="Neuron","5"="Astro","6"="Oligo","7"="Neuron","8"="Oligo","9"="Astro","10"="Astro",
                        "11"="VLMC","12"="Oligo")
  
  fetal.neurons <- subset(fetal,idents="Neuron")
  fetal.neurons <- DietSeurat(fetal.neurons)  
  
  group.by="orig.ident"
  
  fetal.neurons <- NormalizeData(fetal.neurons)
  fetal.neurons<-FindVariableFeatures(fetal.neurons,selection.method = "vst",nfeatures = 6000)
  fetal.neurons <- ScaleData(fetal.neurons)
  fetal.neurons<-RunPCA(fetal.neurons,npcs = 50)
  
  harmony.neurons.reintegrated <- fetal.neurons %>% 
    RunHarmony(group.by, plot_convergence = TRUE)
  
  harmony.neurons.reintegrated <- harmony.neurons.reintegrated %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = .1) %>% 
    identity()
  
  DimPlot(harmony.neurons.reintegrated,label = T,raster = F,pt.size = 1)+scale_color_manual(values = stevens.greenblue())
  
  # Figure 2D
  genes <- c("GAP43","SNAP25","PBX1","NCAM1","SYT1","SOX1","LMX1A","OTX2","EN1","NR4A2","TH","RET","SLC18A2","KCNJ6","CALB1")
  ae.hpsc.genes <- AverageExpression(harmony.neurons.reintegrated,features = genes,group.by = "seurat_clusters")
  
  breaksList = seq(-1, 1, by = .1)
  
  myColor.blue_beige <- colorRampPalette(colors = c( "#969FD9","#BDC2E3", "#F1F2F2","#CE9C5F","#20282F" ))(50)
  pheatmap::pheatmap(ae.hpsc.genes$RNA,scale = "row",cluster_cols = F,cluster_rows = F,color = myColor.blue_beige[25:40],breaks=breaksList)
  
  # Figure 2E
  # ROGUE
  # FIXME: Add score calculation and plot
  
  # Figure 2F
  hpsc.neuron.data$THExpr<-GetAssayData(hpsc.neuron.data,slot = "counts")["TH",]
  hpsc.neuron.data$THBinary<-ifelse(hpsc.neuron.data$THExpr>0,"Pos","Neg")
  
  DimPlot(hpsc.neuron.data,group.by = "THBinary",pt.size = 1.1,label = F,cols = c("grey","brown"))
  
  # Figure 2G
  DimPlot(hpsc.thneuron.data,label=T,group.by = "NewCluster")
  
  # Figure 2H
  Idents(hpsc.thneuron.data) <- hpsc.thneuron.data$NewCluster
  aE <- AverageExpression(hpsc.thneuron.data,features = VariableFeatures(hpsc.thneuron.data)[1:2000])$RNA
  dend <- t(aE) %>% scale %>% dist %>% 
    hclust()
  plot(dend)
  
  # Figure 2I
  ss <- hpsc.thneuron.data@meta.data %>% group_by(Age) %>% sample_n(size=600)
  hh.tmp <- hpsc.thneuron.data[,ss$CellID]
  
  hh.tmp$Age2 <- factor(recode_factor(hh.tmp$Age, "3m"= "6m","6m"="3m"),levels = c("3m","6m","9m","12m"))
  DimPlot(hh.tmp,label = F,raster = F,pt.size = 1.5,split.by = "Age2",ncol = 2)
  
  # Figure 2J
  Idents(hpsc.thneuron.data) <- hpsc.thneuron.data$NewCluster
  fm<-FindAllMarkers(hpsc.thneuron.data,only.pos = T,max.cells.per.ident = 100000,logfc.threshold = .1)
  top.markers <- fm %>% group_by(cluster) %>% top_n(n = 10, wt =  avg_log2FC)
  ds<-hpsc.thneuron.data@meta.data %>% group_by(NewCluster) %>% sample_n(size = 50,replace = T)
  
  hpsc.thneuron.data<-ScaleData(hpsc.thneuron.data,features = top.markers$gene)  
  DoHeatmap(hpsc.thneuron.data[,unique(ds$CellID)],features = unique(top.markers$gene),group.by = "NewCluster")+
    scale_fill_gradientn(colors = c( "#969FD9","#BDC2E3", "#F1F2F2","#CE9C5F","#20282F"))
  
  # Figure 2K
  
  genelist<-read.csv("genelists/genelist_violin_fig2K.txt",header = F)
  
  StackedVlnPlot(hpsc.thneuron.data,features = genelist$V1[1:8])
  ggsave(paste0(outDir,"/vln.1.pdf"),w=4,scale = 2)
  StackedVlnPlot(harmony_human.THneurons,features = genelist$V1[9:18],colors = tableau_color_pal(palette = "Miller Stone")(10) )
  ggsave(paste0(outDir,"/vln.2.pdf"),w=4,scale = 2)
  StackedVlnPlot(harmony_human.THneurons,features = genelist$V1[19:29],colors = tableau_color_pal(palette = "Miller Stone")(10) )
  ggsave(paste0(outDir,"/vln.3.pdf"),w=4,scale = 2)
  
  
  
  
  }
