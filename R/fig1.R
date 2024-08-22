figure1 <- function(outdir,hpsc.file,fetal.file) {
  
  dir.create(outdir)
  dir.create(paste0(outdir,"/figure1/"))
  outdir <- paste0(outdir,"/figure1/")
  
  # Prepare hpsc data
  hpsc.data <- readRDS(hpsc.data)
  hpsc.data$UMAP1<-Embeddings(hpsc.data,reduction = "umap")[,1]
  hpsc.data$UMAP2<-Embeddings(hpsc.data,reduction = "umap")[,2]
  hpsc.data$NamedIdents <- factor(Idents(hpsc.data),
                                                  levels=c("Astrocytes","VLMC","Neurons","Neuron precursors","Oligodendrocytes"))
  
  # Prepare fetal data
  fetal.data<-readRDS(fetal.file)
  fetal.data <- RenameIdents(fetal.data, "0"="Astrocytes","1"="Neuron precursors", "2"="Neurons",
                             "3"="Astrocytes","4"="Neurons","5"="Astrocytes","6"="Oligodendrocytes","7"="Neurons","8"="Oligodendrocytes",
                             "9"="Astrocytes","10"="Astrocytes",
                        "11"="VLMC","12"="Oligodendrocytes")
  fetal.data$NamedIdents<-Idents(fetal.data)
  fetal.data$NamedIdents <- factor(fetal.data$NamedIdents,levels = c("Astrocytes","VLMC","Neuron precursors","Neurons","Oligodendrocytes"))
  
  
  # Figure 1M
  fetal.data$UMAP1<-Embeddings(fetal.data,reduction = "umap")[,1]
  fetal.data$UMAP2<-Embeddings(fetal.data,reduction = "umap")[,2]
  
  ggplot(fetal.data@meta.data[],
         aes(x=UMAP1,y=UMAP2,fill=NamedIdents))+
    geom_point(alpha=.7,pch=21,size=3,colour="#3b3b3b")+theme_cowplot()+scale_fill_manual(values = newpan)&NoLegend()
  
  
  # Figure 1N
  ggplot(hpsc.data@meta.data[sample(1:nrow(hpsc.data@meta.data),15000),],
                aes(x=UMAP1,y=UMAP2,fill=NamedIdents2))+
    geom_point(alpha=.6,pch=21,size=3,colour="#141414")+theme_cowplot()+scale_fill_manual(values = newpan)&NoLegend()
  
  # Figure 1O
  rogue.hes <-  rogue(GetAssayData(hpsc.data),
                      platform = "UMI",
                      labels=hpsc.data$NamedIdents,
                      samples=hpsc.data$orig.ident)
  rogue.boxplot(rogue.hes)
  
  rogue.fetal <-  rogue(GetAssayData(feta.data),
                        platform = "UMI",
                        labels=rep("1",ncol(fetal.data)),
                        samples=fetal.ss$orig.ident)
  rogue.boxplot(rogue.fetal)
  
  # Figure 1P
  f.prop<-data.frame(melt(prop.table(table(fetal.data$NamedIdents))),
                  source="fetal")
  hpsc.prop <- data.frame(melt(prop.table(table(hpsc.data$NamedIdents))),
                          source="hpsc")
  
  rbind(f.prop,hpsc.prop) %>% ggplot(aes(x=source,y=value*100,fill=Var.1))+
    geom_bar(stat="identity")+ylab("cells (%)")+xlab("")+
    theme_cowplot()+scale_fill_manual(values = newpan)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))&NoLegend()
  
  # Figure 1Q
  hpsc.data<-FindVariableFeatures(hpsc.data,nfeatures = 12000)
  hpsc.data<-ScaleData(hpsc.data)
  hpsc.data<-RunPCA(hpsc.data)
  anchors <- FindTransferAnchors(reference = fetal.data, query = hpsc.data, 
                                          dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = fetal.data$NamedIdents, 
                              dims = 1:30)
  hpsc.data <- AddMetaData(hpsc.data, metadata = predictions)
  
  colfunc <- colorRampPalette(c("#F3F5F9","#919AA7", "#D5A05A"))
  FeaturePlot(hpsc.data,features = "prediction.score.max",min.cutoff = "q2")+
    scale_colour_gradientn(colours = colfunc(20))+ggtitle("")
  
  # Figure 1R
  VlnPlot(hpsc.data,"prediction.score.max",
          pt.size = 0)&scale_fill_manual(values = newpan)&NoLegend()
  
  # Figure 1S
  
  pp<-prop.table(table(hpsc.data$NamedIdents,hpsc.data$predicted.id),1)
  
  ppm <- melt(pp)
  
  colnames(ppm)[3]<-"Percentage"
  ppm$Percentage<-ppm$Percentage*100
  
  ggplot(ppm,aes(x=Var.1,y=Var.2,size=Percentage,col=Var.1))+geom_point()+
    theme_cowplot()+scale_size(range=c(0,30))+xlab("hES")+ylab("Fetal")+
    scale_color_manual(values = newpan) + guides(col = FALSE)
  
  
  
  }
