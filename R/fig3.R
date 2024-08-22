figure3 <- function(outdir,hpsc.thneuron.file="processed_data/hpsc.human_harmony.THneurons.rds"){
  
  dir.create(outdir)
  dir.create(paste0(outdir,"/figure3/"))
  outdir <- paste0(outdir,"/figure3/")
  
  # Prepare hpsc data
  hpsc.thneuron.data <- readRDS(hpsc.thneuron.file)
  hpsc.thneuron.data$NewCluster <- factor(paste0("DA-",as.integer(hpsc.thneuron.data$seurat_clusters)))
  hpsc.thneuron.data$cellbc_sample <- colnames(hpsc.thneuron.data)
  
  # prepare mcherry count data
  mch <- read.csv("processed_data/mcherry_counts.csv")
  
  mch.full <- merge(hpsc.thneuron.data@meta.data,mch,by="cellbc_sample",all.x=T)
  rownames(mch.full)<-mch.full$cellbc_sample
  
  mch.full[is.na(mch.full$mCherry_counts),"mCherry_counts"]<-0
  hpsc.thneuron.data<-AddMetaData(hpsc.thneuron.data,mch.full)
  
  
  # Figure 3G
  Idents(hpsc.thneuron.data)<-hpsc.thneuron.data$tracing
  fm.volc<-FindMarkers(hpsc.thneuron.data,min.pct = 0,logfc.threshold = 0,ident.1="Striatum",ident.2="PFC")
  genes.volc <- read.csv("genelists/3g.volcano.txt",header = F)

  keyvals.shape <- ifelse(
    rownames(fm.volc) %in% genes.volc$V1[1:9], "darkgreen",
    ifelse(rownames(fm.volc) %in% genes.volc$V1[10:19], "darkblue",
           "grey"))
  
  names(keyvals.shape)<- NA
  names(keyvals.shape)[keyvals.shape == "darkgreen"] <- 'Striatum'
  names(keyvals.shape)[keyvals.shape == "darkblue"] <- 'PFC'
  
  EnhancedVolcano(fm.volc,
                  lab = rownames(fm.volc),
                  x = 'avg_log2FC',selectLab = genes.volc$V1,
                  y = 'p_val',pCutoff = 1,FCcutoff = .5,ylim = c(0,8),xlim = c(-2,2),subtitle = "",
                  colCustom = keyvals.shape,drawConnectors = TRUE,title = "",
                  pointSize = ifelse(rownames(fm.volc)%in% genes.volc$V1,8,1))
  
  # Figure 3H
  eg <- bitr(rownames(fm.volc[fm.volc$avg_log2FC > 0 & fm.volc$p_val<0.05,]),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  univ <- bitr(rownames(hpsc.thneuron.data),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  
  ego <- enrichGO(gene          = eg$ENTREZID,
                  OrgDb         = org.Hs.eg.db,universe = univ$ENTREZID,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.1,
                  readable      = TRUE)
  edox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
  cnetplot(edox, categorySize="pvalue")
  
  # Figure 3I
  eg <- bitr(rownames(fm.volc[fm.volc$avg_log2FC < 0 & fm.volc$p_val<0.05,]),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  univ <- bitr(rownames(hpsc.thneuron.data),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  
  ego <- enrichGO(gene          = eg$ENTREZID,
                  OrgDb         = org.Hs.eg.db,universe = univ$ENTREZID,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.1,
                  readable      = TRUE)
  edox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
  cnetplot(edox, categorySize="pvalue")
  
  
  }
