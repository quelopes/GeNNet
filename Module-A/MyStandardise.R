myStandardise = function(data){
  
  # data  = mat.exp
  
  if(ncol(data) > 2 ){
    for (i in 1:nrow(data)) {
      data[i, ] <- (data[i, ] - mean(data[i, ], na.rm = TRUE))/sd(data[i,], na.rm = TRUE)
    }
    data
  }
  if(ncol(data) == 2){
    data
    
  }
}
