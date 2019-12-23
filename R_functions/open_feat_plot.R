

feature_open <- function(dir = "R:/Akina/Raw SC data/human_1702/59d_SI-3A-E3/outs/filtered_gene_bc_matrices/GRCh38", obj = NULL, GOI = c("VSX1", "GNAO1", "GRM6", "TRPM1"), title = "D59", th = 0){
  
  if(is.null(obj)){
    Se <- Read10X(data.dir = dir)
    Se <- CreateSeuratObject(counts =  Se)
    Se <- NormalizeData(Se, verbose=T)
    Se <- ScaleData(Se, verbose = T)
  }
  else
    {Se <- obj}
  print(Se)
  Se <- FindVariableFeatures(Se, selection.method="vst", nfeatures=2000)
  Se <- RunPCA(Se, features = VariableFeatures(object = Se), verbose = T, npcs = 70)
  Se <- FindNeighbors(Se)
  Se <- FindClusters(Se)
  Se <- RunUMAP(Se, dims = 1:70)
  p = FeaturePlot(Se, features = GOI, min.cutoff =  th, order = T) + ggtitle(title)
  print(p)
  return(Se)
}