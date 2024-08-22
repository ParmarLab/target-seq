figure3 <- function(outdir,hpsc.neuron.file){
  
  dir.create(outdir)
  dir.create(paste0(outdir,"/figure3/"))
  outdir <- paste0(outdir,"/figure3/")
  
  # Prepare hpsc data
  hpsc.thneuron.data <- readRDS("processed_data/hpsc.human_harmony.THneurons.rds")
  hpsc.thneuron.data$NewCluster <- factor(paste0("DA-",as.integer(hpsc.thneuron.data$seurat_clusters)))
  hpsc.thneuron.data$cellbc_sample <- colnames(hpsc.thneuron.data)
  
  # prepare mcherry count data
  mch <- read.csv("processed_data/mcherry_counts.csv")
  
  mch.full <- merge(hpsc.thneuron.data@meta.data,mch,by="cellbc_sample",all.x=T)
  rownames(mch.full)<-mch.full$cellbc_sample
  
  mch.full[is.na(mch.full$mCherry_counts),"mCherry_counts"]<-0
  hpsc.thneuron.data<-AddMetaData(hpsc.thneuron.data,mch.full)
  hpsc.thneuron.data$mCherry_counts_binary<- ifelse(hpsc.thneuron.data$mCherry_counts > 0,1,0)
  
  # Figure 4A
  pfc<-Nebulosa::plot_density(hpsc.thneuron.data[,hpsc.thneuron.data$tracing != "Striatum"],features = "mCherry_counts_binary")
  striatum<-Nebulosa::plot_density(hpsc.thneuron.data[,hpsc.thneuron.data$tracing != "PFC"],features = "mCherry_counts_binary")
  pfc+striatum
  
  # Figure 4B
  prop_test <- sc_utils(hpsc.thneuron.data)
  
  prop_test <- permutation_test(
    prop_test, cluster_identity = "NewCluster",
    sample_1 = "PFC", sample_2 = "Striatum",
    sample_identity = "tracing"
  )
  permutation_plot(prop_test)+theme(legend.position = "top")
  
  # Figure 4C
  Idents(hpsc.thneuron.data)<-hpsc.thneuron.data$NewCluster
  StackedVlnPlot(hpsc.thneuron.data[,hpsc.thneuron.data$NewCluster %in% c("DA-7","DA-2")],
                 features = c("KCNJ6","SLC17A6","SNCA","SLC18A2","OTX2","FOXP2","LRP4","SEZ6L","LMO3","TCF4")
                 ,colors = c("darkblue","darkgreen"))
  
  # Prepare Figure 4D-I
  
  markers_pfc_str <- wilcoxauc(hpsc.thneuron.data,"tracing", assay = "data", 
                                      groups_use = c('Striatum', 'PFC' ))
  gbm <- getBM(attributes = c('hgnc_symbol', 'entrezgene_id'), 
               filters = 'hgnc_symbol', 
               values = markers_pfc_str$feature, 
               mart = ensembl)
  
  markers_pfc_str <- merge(markers_pfc_str,gbm,by.x="feature",by.y="hgnc_symbol")
  
  markers_pfc_str <- markers_pfc_str %>% group_by(entrezgene_id) %>% top_n(n=1,wt = -logFC)
  markers_pfc_str <- markers_pfc_str[!duplicated(markers_pfc_str$entrezgene_id),]
  
  fcs <- markers_pfc_str$logFC
  names(fcs)<-as.character(markers_pfc_str$entrezgene_id)
  fcs<-fcs[!duplicated(names(fcs))]
  
  m_df<-  msigdbr(species = "Homo sapiens")
  
  m_list_c2<- m_df %>% filter(gs_cat == "C2", gs_subcat=="CP:REACTOME") %>%
    split(x = .$gene_symbol, f= .$gs_name)
  
  ranks_PFC<- markers_pfc_str %>%
    filter(group == "PFC") %>%
    arrange(desc(auc)) %>%
    dplyr::select(feature, auc) %>% tibble::deframe()
  
  ranks_Str<- markers_pfc_str %>%
    filter(group == "Striatum") %>%
    arrange(desc(auc)) %>%
    dplyr::select(feature, auc) %>%
    tibble::deframe()
  
  fgsea_resPFC<-fgsea(pathways = m_list_c2, stats = ranks_PFC)
  fgsea_resStr<-fgsea(pathways = m_list_c2, stats = ranks_Str)

  # Figure 4D
  ReactomePA::viewPathway("MAPK family signaling cascades", 
               readable = TRUE, 
               foldChange = fcs)
  # Figure 4E
  plotEnrichment(m_list_c2[["REACTOME_MAPK_FAMILY_SIGNALING_CASCADES"]],
                 ranks_Str) + labs(title="REACTOME_MAPK_FAMILY_SIGNALING_CASCADES")
  
  # Figure 4F
  genes.volc.signi <-as.vector(unlist(fgsea_resStr[fgsea_resStr$pathway =="REACTOME_MAPK_FAMILY_SIGNALING_CASCADES","leadingEdge"]))

  ds<-hpsc.thneuron.data@meta.data %>% group_by(tracing) %>% sample_n(size = 35,replace = F)
  ds <- ds[ds$tracing!="Not traced",]
  
  hpsc.thneuron.data <- ScaleData(hpsc.thneuron.data,features = genes.volc.signi)
  Idents(hpsc.thneuron.data)<-hpsc.thneuron.data$tracing
  dhm<-DoHeatmap(hpsc.thneuron.data[,unique(ds$CellID)],
                 features = genes.volc.signi)
  dhm + scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
  
  # Figure 4G
  viewPathway2("Signaling by WNT", 
               readable = TRUE, 
               foldChange = fcs)&scale_fill_viridis()
  
  # Figure 4H
  plotEnrichment(m_list_c2[["REACTOME_SIGNALING_BY_WNT"]],
                 ranks_Str) + labs(title="REACTOME_SIGNALING_BY_WNT")
  
  # Figure 4I
  genes.volc.signi.wnt <-as.vector(unlist(fgsea_resPFC[fgsea_resPFC$pathway =="REACTOME_SIGNALING_BY_WNT","leadingEdge"]))
  
  hpsc.thneuron.data <- ScaleData(hpsc.thneuron.data,features = genes.volc.signi.wnt)
  dhm<-DoHeatmap(hpsc.thneuron.data[,unique(ds$CellID)],
                 features = genes.volc.signi.wnt)
  dhm + scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
  
  # Figure 4J
  
  
  
  }
