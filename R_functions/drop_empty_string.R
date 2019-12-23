drop_empty <- function(x){
  new_x = NULL
  for(entry in x){
    if (entry != ""){
      new_x = c(new_x, entry)
    }
  }
  return(new_x)
}