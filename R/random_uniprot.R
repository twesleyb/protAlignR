# Function to get some random proteins.
random_uniprot <- function(n){
	library(org.Mm.eg.db)
	db <- org.Mm.egUNIPROT
	all_uniprot <- unlist(as.list(db[mappedkeys(db)]))
	rand_prot <- sample(all_uniprot,n)
	return(rand_prot)
}
