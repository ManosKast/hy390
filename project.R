
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


calculate_pairwise_differences <- function(genome_list) {
  n <- length(genome_list)
  dx <- 1
  differences_list <- numeric(n * (n - 1) / 2)
  
  genome_list <- lapply(genome_list, function(g) strsplit(g, "")[[1]])
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Calculate differences using pre-split genomes
      differences_list[dx] <- sum(genome_list[[i]] != genome_list[[j]])
      dx <- dx + 1
    }
  }
  
  return(differences_list)
}



calculate_average_number_of_pairwise_differences <- function(differences_list){
  total_genomes <- length(differences_list)
  k<-sum(differences_list)/choose(total_genomes, 2)
  return(k)
}


return_500_largest_k <- function(k){
  tmp <- sort(k, decreasing=TRUE)
  top_500_k <- tmp[1:500]
  return(top_500_k)
}

index <- 1
genome_lists <- extract_genome_lists("ms_sim_final.out")
k <- numeric(length(genome_lists))
pairwise_differences <- numeric()

k_500 <- numeric(500)
time_start <- Sys.time()
for(i in seq_along(genome_lists)){
  pairwise_differences <- calculate_pairwise_differences(genome_lists[[i]])
  k[i] <- calculate_average_number_of_pairwise_differences(pairwise_differences)
  #print(index)
  #index <- index + 1
}

k_500 <- return_500_largest_k(k)
print(k_500)

time_end <- Sys.time()

print(time_end - time_start)

