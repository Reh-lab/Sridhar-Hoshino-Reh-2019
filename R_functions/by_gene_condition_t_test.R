by_gene_t.test <- function(SerautObj, genes.use, group.by=NULL, assay = "RNA", comp = NULL, alpha_start = .05, Bonferroni = T) {
  ### Initialize values 
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  
  if (length(comp) != 2){
    warning("This function only compares two treatments at this time. Returning -1.")
    return(-1)
  }
  
  # a list of all possible splits
  #groups = unlist(as.list(unique(rgcsubset[["time"]])))
  

  #first we want the gene vector for each of our genes
  for (gene in genes.use){
    
    
    #then we constuct the gene vectors
    i_cells = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
    i_vector = SerautObj@assays[[assay]]@data[gene, i_cells]
    
    
    j_cells = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
    j_vector = SerautObj@assays[[assay]]@data[gene, j_cells]
    
    
    
    #preform t-test
    t_out = t.test(i_vector, j_vector)

    # assign condition code
    cond = paste(comp[1], comp[2], sep = "_")
    
    
    #constuct outs
    
    condition.out <- c(condition.out, cond)
    stat.out <- c(stat.out, t_out[["statistic"]])
    p_val.out <- c(p_val.out, t_out[["p.value"]])
    gene.out <- c(gene.out, gene)
    
    
  }
  
  if (Bonferroni == T){
    #  Bonferroni correction
    new_alpha = alpha_start/(2*length(genes.use))
    
    # report alpha
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    
    # create significance vector
    sig_out = p_val.out < new_alpha
    
    #create df
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)
  }
  
  else{
    #create df
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out)
  }
  
  
  return(dfOUT)
}