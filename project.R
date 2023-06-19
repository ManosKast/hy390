#==============================================================================

#                 a, b, c and e related functions


calculate_a1 <- function(number) {
  a1 <- 0
  
  for(i in 1:number){
    a1 <- a1 + 1/i
  }
  
  return(a1)
}


calculate_a2 <- function(number) {
  a2 <- 0
  
  for(i in 1:number){
    a2 <- a2 + 1/i^2
  }
  
  return(a2)
}


calculate_b1 <- function(n) {
  b1 <- (n+1)/(3*(n-1))
  return(b1)
}


calculate_b2 <- function(n) {
  b2 <- 2*(n^2 + n + 3) / (9*(n * (n - 1)))
  return(b2)
}



calculate_c1 <- function(a1, b1) {
  c1 <- b1 - 1/a1
  return(c1)
}



calculate_c2 <- function(a1, a2, b2, n) {
  c2 <- b2 - (n+2)/(a1 * n) + a2/a2^2
  return(c2)
}


calculate_e1 <- function(a1, c1) {
  e1 <- c1/a1
  return(e1)
}


calculate_e2 <- function(a1, a2, c2) {
  e2 <- c2 / (a1^2 + a2)
  return(e2)
}


#==============================================================================


#               GET GENOME LIST


extract_genome_lists <- function(filename) {
  genomes <- readLines(filename)
  genome_collections <- list()
  genome_strings <- character()
  
  for (genome in genomes) {
    if(genome != ""){
      genome_strings <- c(genome_strings, genome)
    }
    else {
      genome_collections <- c(genome_collections, list(genome_strings))
      genome_strings <- character()
    }
  }  
  
  return(genome_collections)
}

extract_final_list <- function(filename) {
  genomes <- readLines(filename)
  genome_collections <- list()
  genome_strings <- character()
  
  for (genome in genomes) {
    genome_strings <- c(genome_strings, genome)
  }  
  
  genome_collections <- c(genome_collections, list(genome_strings))

  return(genome_collections)
}


#==============================================================================

#                      K-RELATED FUNCTONS

calculate_pairwise_differences <- function(genome_list) {
  n <- length(genome_list)
  dx <- 1
  differences_list <- numeric(n * (n - 1) / 2)
  genome_list <- lapply(genome_list, function(g) strsplit(g, "")[[1]])

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      differences_list[dx] <- sum(genome_list[[i]] != genome_list[[j]])
      dx <- dx + 1
    }
  }
  
  return(differences_list)
}



calculate_k <- function(genome_list){
  differences_list <- calculate_pairwise_differences(genome_list)
  total_genomes <- length(differences_list)
  k<-sum(differences_list)/total_genomes
  return(k)
}


#==============================================================================

#                 W-RELATED FUNCTIONS

calculate_w <- function(genome_list) {
  w <- nchar(genome_list[[1]]) / calculate_a1(length(genome_list)-1)
  return(w)
}



#==============================================================================

#               D-RELATED FUNCTIONS

calculate_D <- function(genome_list, k, w) {
  n <- length(genome_list)
  S <- nchar(genome_list[[1]])
  
  a1 <- calculate_a1(n-1)
  a2 <- calculate_a2(n-2)
  
  b1 <- calculate_b1(n)
  b2 <- calculate_b2(n)
  
  c1 <- calculate_c1(a1, b1)
  c2 <- calculate_c2(a1, a2, b2, n)
  
  e1 <- calculate_e1(a1, c1)
  e2 <- calculate_e2(a1, a2, c2)
  
  expr <- e1 * S + e2 * S * (S - 1)
  D <- (k - w) / sqrt(expr)
  
  return(D)
}





return_500_largest_values <- function(k){
  tmp <- sort(k, decreasing=TRUE)
  top_500_k <- tmp[1:500]
  return(top_500_k)
}




#       ERWTHMATA 1 KAI 2
#MENEI TO D0, w0, k0.

final_genome_list <- extract_final_list("ms_obs_final.out")
k0 <- calculate_k(final_genome_list[[1]])
w0 <- calculate_w(final_genome_list[[1]])
D0 <- calculate_D(final_genome_list[[1]], k0, w0)


genome_lists <- extract_genome_lists("ms_sim_final.out")
k <- numeric(length(genome_lists))
w <- numeric(length(genome_lists))
D <- numeric(length(genome_lists))


time_start <- Sys.time()
for(i in seq_along(genome_lists)){
  k[i] <- calculate_k(genome_lists[[i]])
  
  w[i] <- calculate_w(genome_lists[[i]])
  
  D[i] <- calculate_D(genome_lists[[i]], k[i], w[i])
}



time_end <- Sys.time()

print(time_end - time_start)



#     ERWTHMA 3

standarise_vector <- function(vector) {
  i <- 1
  standarised_vector <- numeric(length(vector))
  
  average_value <- mean(vector)
  std_deviation <- sd(vector)
  
  for(value in vector){
    standarised_vector[i] <- (value - average_value)/std_deviation
    i <- i + 1
  }
  
  return(standarised_vector)
}


standarise_obs_vector <- function(obs_value, sim_vector){

  average_value <- mean(sim_vector)
  std_deviation <- sd(sim_vector)

  standarised_vector <- (obs_value - average_value)/std_deviation

  
  return(standarised_vector)
}


standarised_k <- standarise_vector(k)
standarised_w <- standarise_vector(w)
standarised_D <- standarise_vector(D)


standarised_k0 <- standarise_obs_vector(k0, k)
standarised_w0 <- standarise_obs_vector(w0, w)
standarised_D0 <- standarise_obs_vector(D0, D)


#     ERWTHMA 4

calculate_Euclidean_distance <- function(k0, w0, D0, k, w, D){
  len <- length(k)
  d <- numeric(len)
  for(i in 1:len) {
    expr <- (D0 - D[[i]])^2 + (w0 - w[[i]])^2 + (k0 - k[[i]])^2
    d[i] <- expr
  }
  
  d <- sqrt(d)
  return(d)
}


d <- calculate_Euclidean_distance(standarised_k0, standarised_w0, standarised_D0,
                                  standarised_k, standarised_w, standarised_D)





