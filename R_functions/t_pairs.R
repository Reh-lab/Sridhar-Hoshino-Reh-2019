t_pairs <- function(pairs_trt, df, genes, v = F, alpha = .05){
  
  t_df = NULL
  for (i in 1:length(pairs_trt)){
    print(paste("Pair number", i))
    print(pairs_trt[[i]][[1]])
    print(pairs_trt[[i]][[2]])
    tres = by_gene_t.test(df, group.by="time_condense_type",
                          genes.use = genes, comp = c(pairs_trt[[i]][[1]], pairs_trt[[i]][[2]]),
                          Bonferroni = F)
    if (v==T){
      print(tres)
    }
    t_df = rbind(tres, t_df)
  }
  
  # do Bonferroni for merged data
  n_alpha = alpha / length(t_df$p_val)
  #make significance
  t_df$sig = t_df$p_val < n_alpha
  return(t_df)
}