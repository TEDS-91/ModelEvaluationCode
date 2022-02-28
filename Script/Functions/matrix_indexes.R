matriz_indexes <- function(matriz){
  
  n <- nrow(matriz)
  
  ix <- matrix(NA, n * n, 2)
  
  k <- 0
  
  for (i in (1:n)){
    for (j in (1: n)){
      k <- k + 1 
      ix[k, 1]<- i
      ix[k, 2]<- j
    }
  }
  
  o <- rep(NA, nrow(ix))
  
  o <- matriz[ix]
  
  out <- cbind(ix, o)
  
  colnames(out) <- c("RowInd", "ColInd", "Values")
  
  list(matriz = out,
       funcao = "Funcao criada por Talita Santana.")
  
}