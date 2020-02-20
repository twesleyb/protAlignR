#!/usr/bin/env Rscript

# Alignment of protein sequences with EMBOSS::water algorithm.

library(getPPIs)
library(seqinr)


# Align Fam81a to proteins of interest.
protA <- getPPIs::mapIds("Fam81a",from="symbol",to="uniprot",species="mouse")

# Get uniprot ids of proteins of interest.
prots <- c("Cacng8","Kcnq4","Shisa7","Ctps","Prrt1","Strip1")
uniprot <- getPPIs::mapIds(prots,from="symbol",
			   to="uniprot",species="mouse")
uniprot[1] <- "Q8VHW2"

# Loop to align proteins.
result <- list()
for (prot in uniprot) {
       result[[prot]] <- algnWater(protA,prot)
}
names(result) <- paste(names(uniprot),uniprot,sep="|")

res <- result[["Kcnq4|Q9JK97"]]

# Exctract scores.
scores <- sapply(result,function(x) x$Summary$Score)

# Function to get some random proteins.
random_uniprot <- function(n){
	library(org.Mm.eg.db)
	db <- org.Mm.egUNIPROT
	all_uniprot <- unlist(as.list(db[mappedkeys(db)]))
	rand_prot <- sample(all_uniprot,n)
	return(rand_prot)
}

# Number of permutations required for significance.
alpha <- 0.05/length(uniprot)
nperm <- NetRep::requiredPerms(alpha, alternative = "greater")

# Loop to align Fam81a with random proteins.
rand_prots <- random_uniprot(nperm)
rand_result <- list()
pbar <- txtProgressBar(max=length(rand_prots),style=3)
for (i in seq_along(rand_prots)) {
	setTxtProgressBar(pbar,i)
	prot <- rand_prots[i]
	rand_result[[i]] <- algnWater(protA,prot)
}
close(pbar) ; message("\n")

# Examine the result.
length(rand_result)
rand_scores <- sapply(rand_result,function(x) x$Summary$Score)

# Compare observed scores to random scores.
p <- sapply(scores, function(x) sum(rand_scores>x)/length(rand_scores))
q <- p*nperm

# Function to perform water alignment.
algnWater <- function(protA,protB){
	suppressPackageStartupMessages({
		library(seqinr)
	})
	# Execute run-water to align proteins.
	cmd <- paste("./run-water",protA,protB)
	result <- system(cmd,intern=TRUE)
	# Parse the result.
	delim <- "#======================================="
	data <- splitAt(result,grep(delim,result))
	names(data) <- c("parameters","summary","alignment")
	# Clean up the result summary.
	raw <- data[[2]]
	temp <- strsplit(raw,":")
	temp <- temp[which(sapply(temp,length)==2)]
	params <- vector("list",length(temp))
	names(params) <- gsub("# ","",sapply(temp,"[",1))
	names(params)[2] <- "Protein_1"
	names(params)[3] <- "Protein_2"
	vals <- trimws(sapply(temp,"[",2))
	params[] <- vals
	idx <- c(1,5,6,7,11)
	params[idx] <- lapply(params[idx],as.numeric)
	# Collect the protein sequences.
	fasta_files <- paste0(c(protA,protB),".fasta")
	#names(fasta_files) <- c(protA,protB)
	names(fasta_files) <- fasta_files
	prot_seq <- sapply(fasta_files, function(fasta) {
			   read.fasta(fasta,seqtype="AA", 
				      as.string=TRUE,seqonly=TRUE) })
	unlink(fasta_files)
	# Clean up the alignment.
	raw <- data[[3]]
	max_nchar <- 13
	a <- substr(params$Protein_1,1,max_nchar)
	b <- substr(params$Protein_2,1,max_nchar)
	idx <- mapply(zip,grep(a,raw),grep(b,raw))
	idx <- unlist(lapply(idx,function(x) c(x[1],x[1]+1,x[2])))
	algn <- raw[idx]
	# Collect results.
	results <- list()
	results[["Input"]] <- prot_seq
	results[["Summary"]] <- params
	results[["Alignment"]] <- algn
	results[["stdout"]] <- result
	# Return the results.
	return(results)
}
