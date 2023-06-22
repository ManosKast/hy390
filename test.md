Το παρόν σκριπτ αποσκοπεί να αναλύσει δεδομένα γονιδίων μέσω της μεθόδου του Tajima.
Αρχικά, ορίζονται αρκετές βοηθητικές συναρτήσεις για τον υπολογισμό της μεθόδου του Tajima D και της μεθόδου του Watterson --calculate_a1, calculate_a2, calculate_b1, calculate_b2, calculate_c1, calculate_c2, calculate_e1, calculate_e2.

Έπειτα, ορίζονται δύο συνρτήσεις που διαβάζουν και αναλύουν τα δοσμένα γονίδια, δημιουργόντας λίστες γονιδίων από τη προσομοίωση και τη παρατήρηση --extract_genome_lists, extract_final_list.

Οι βασικότερες συναρτήσεις του script είναι οι calculate_k, calculate_w και calculate_D, οι οποίες υπολογίζουν τον μέσο όρο των διαφορών ανά ζεύγος, τη μέθοδο του Watterson και τη μέθοδο του Tajima D, αντίστοιχα.

Ξεκινώντας το πρόγραμμα, αποθηκεύω στις μεταβλητές k0, w0 και D0 τις τιμές που παράγονται συναρτήσει των παρατηρούμενων γονίδιων, που βρίσκονται εντός του αρχείου ms_obs_final.out.
Αντίστοιχα, αποθηκεύω στις μεταβλητές k, w, D τις τιμές που παράγονται συναρτήσει των γονιδίων της προσομοίωσης που βρίσκονται εντός του αρχείου ms_sim_final.out.

Έπειτα, για να μπορέσουμε να συγκρίνουμε δικαίως τις τιμές των αποτελεσμάτων των παραπάνω μεταβλητών, τυποποιώ μέσω της standarise_vector τις τιμές των αποτελεσμάτων της προσομοίωσης και μέσω της standarise_obs_vector τις τιμές των αποτελεσμάτων των μεταβλητών των k0, w0 και D0.
Χρησιμοποιώντας αυτές τις τυποποιημένες τιμές, μέσω του πυθαγορείου θεωρήματος υπολογίζω τη διαφορά 'αποστάσεως'  μεταξύ τ ων προσομοιωμένων τιμώ ν και των παρατηρούμενων --συνάρτηση calculate_Euclidean_distance.

Προχωρώντας, μέσω της συνάρτησης get_500_smallest_values_indexes, λαμβάνω τους δείκτες με τις 500 μικρότερες τιμές, καθώς δηλώνουν τη καλύτερη δυνατή προσέγγιση.
Έπειτα, μέσω των 500 αυτών δεικτών επιλέγω τις 500 παραμέτρους που βρίσκονται εντός του αρχείου pars_final.txt --συνάρτηση get_parameter_values.
Αναζητώ τη μέση τιμή και τη διάμεσο των 500 αυτών παραμέτρων και τέλος παράγω ένα ιστόγραμμα και ένα διάγραμμα πυκνότητας για τις τιμές αυτές.

![R_project_graph](https://github.com/ManosKast/hy390/assets/92722366/c2e95b7d-87b7-47ad-9cc1-7de528b292de)
![R_project_mean_median](https://github.com/ManosKast/hy390/assets/92722366/288b6f88-48cf-4e83-94ad-0b27db3ac300)

Παρατηρώ ότι η μέση τιμή των παραμέτρων είναι 109.5133 και η διάμεσος των παραμέτρων είναι 106.6795.
Παρατηρώντας το διάγραμμα πυκνότητας, παρατηρώ ότι το πλήθος των παραμέτρων σχεδόν παράγει μία δεξιοστροφη κανονική κατανομή.

Η τιμή διαμέσου με τη μέση τιμή είναι σχεδόν ισοδύναμες, συνεπώς πράγματι το διάγραμμα πυκνότητας θυμίζει ένα κανονικό διάγραμμα με κέντρο το τη μέση τιμή και αφού η μέση τιμή είναι μεγαλύτερη από τη τιμή διαμέσου η λόξα αυτή επιβεβαιώνεται.
Από τη προσομοίωση και το διάγραμμα πυκνότητας συνάγω ότι το διάγραμμα πλησιάζει ένα δεξιόστροφο διάγραμμα κανονικής κατανομής με κέντρο τη μέση τιμή. Η παρατήρηση αυτή μπορεί να επιβεβαιωθεί και από το γεγονός ότι η μέση τιμή είναι μεγαλύτερη από τη τιμή διάμεση τιμή, αλλά παρ'όλα αυτά είναι σχεδόν ισοδύναμες.
Αυτό υπονοεί πως ο μέσος ετήσιος ρυθμός ανάπτυξης θα ήταν ασφαλές να υποθέσουμε ότι θα είναι κοντά στις 106-110 μονάδες, βάσει των συνθηκών των προσομοιώσεων μας.
Η δεξιά λόξα επιδυκνύει πως υπάρχει ένας αριθμός ενδεχομένων ο COVID να αναπτυχθεί με αρκετά ταχύτερους ρυθμούς, αλλά αν και πιθανό ενδεχόμενο είναι λιγότερο πιθανό κάτι τέτοιο να συμβεί.

Κώδικας R:

```R
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


for(i in seq_along(genome_lists)){
  k[i] <- calculate_k(genome_lists[[i]])
  
  w[i] <- calculate_w(genome_lists[[i]])
  
  D[i] <- calculate_D(genome_lists[[i]], k[i], w[i])
}




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




#             ERWTHMA 5


get_500_smallest_values_indexes <- function(vector){
  sorted_indexes <- order(vector)
  smallest_indexes <- sorted_indexes[1:500]
  return(smallest_indexes)
}

smallest_500_indexes <- get_500_smallest_values_indexes(d)


#             ERWTHMA 6


# indexes variable contains the indexes that correspond to the value's line.
# Given those indexes, append to a vector these values and return it.
get_parameter_values <- function(filename, indexes){
  lines <- readLines(filename)
  
  values <- as.numeric(lines[indexes])
  return(values)
}


parameters <- get_parameter_values("pars_final.txt", smallest_500_indexes)


#         ERWTHMA 7

mean_parameters <- mean(parameters)
median_parameters <- median(parameters)

cat("Mean: ")
print(mean_parameters[[1]])
cat("Median: ")
print(median_parameters[[1]])


#         ERWTHMA 8

dens <- density(parameters)

par(mfrow = c(1, 2))
hist(parameters, breaks = 10, col = "skyblue", border = "black", xlab = "Parameters", ylab = "Frequency", main = "Histogram of Parameters")
plot(dens, main = "Density Plot of Parameters")
```
