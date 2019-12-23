recolor_plot <- function(Se, csv = "Z:/People Folders/Connor/From Akshaya/recolor.csv", show_orig_plot = F, start_msg = T, other_name = "Other", label = F, ax = F, pt.size = 1){
  
  if (start_msg){
    print("Before starting this function there are a few things you should know. this function takes a path to a csv and a seurat object. It assumes that whatever you want colored on the dimplot is in the idents slot. It also assumes that the csv contains two or more columns, where the first column contains the ident names, and the second column contains color codes or names (ex. grey40 or #AEAEE5). It will tell you which idents it does not find in your csv, you can add these to the csv or leave them alone. The values it cannot find will be assigned the identy Other which must be in the csv or this will break. I can be more strict in my checks so the function is harder to break if you would like. Turn this message off by adding the argument , start_msg = F")
  }
  
  Se_tmp = Se
  
  library(dplyr)
  library(ggplot2)
  
  recolor = read.csv(csv)
  
  # get levels
  levels(Idents(Se_tmp))
  
  #save idents
  Se_tmp$type <- Idents(Se_tmp)
  Se_tmp$type2 <- as.character(Se_tmp$type)
  
  ill_names = c()
  #pill_names = c() this could be partial matches
  
  
  for(i in as.character(unique(Se_tmp$type))){
    if(!(i %in% recolor[,1])){
      message(paste("The following identity was not found in the color map", i))
      ill_names = c(ill_names, i)
      message("Reclassifying as other")
      message(" ")
      Se_tmp@meta.data[Se_tmp@meta.data$type == i, "type2"] <- other_name
    }
  }
  Idents(Se_tmp) <- "type2"
  
  recolor = recolor[recolor[,1] %in% levels(Idents(Se_tmp)),]
  row.names(recolor) = recolor[,1]
  recolor = recolor[levels(Idents(Se_tmp)),]
  #recolor3
  
  
  #message(paste0("The following names will be reclassified as other ", ill_names)) prints funny
  #pill_names
  
  #as.character(recolor[,1])
  #as.character(recolor3[,1])
  #levels(Idents(Se_tmp))
  print(DimPlot(Se_tmp, label = label, pt.size = pt.size, cols = as.character(recolor[,2])) + NoAxes())
  
  if(show_orig_plot){
    print(DimPlot(Se_tmp, group.by = "type", label = label) + NoAxes())
  }
  
}
