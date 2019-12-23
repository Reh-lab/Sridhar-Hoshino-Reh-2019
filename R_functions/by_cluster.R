
by_cluster <- function(){
  num_cells = length(colnames(fetal_recolored))
  by_type = table(fetal_recolored@active.ident)
  
  n = NULL
  fracT = NULL
  for (t in names(by_type)){
    n = c(t, n)
    val = by_type[[t]]/num_cells
    fracT = c(val, fracT)
  }
  names(fracT) = n
  data.frame("Total" = fracT)
}