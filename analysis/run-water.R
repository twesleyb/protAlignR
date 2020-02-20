#!/usr/bin/env Rscript

# Alignment of protein sequences with EMBOSS::water algorithm.

library(getPPIs)
library(seqinr)

devtools::load_all()

here <- getwd()
root <- dirname(here)
datadir <- file.path(root,"data")

# EBLOSUM62
dm <- readRDS(file.path(datadir,"EBLOSUM62.RData"))

# Align Fam81a to proteins of interest.
protA <- getPPIs::mapIds("Fam81a",from="symbol",to="uniprot",species="mouse")

# Get uniprot ids of proteins of interest.
prots <- c("Cacng8","Kcnq4","Shisa7","Ctps","Prrt1","Strip1")
uniprot <- getPPIs::mapIds(prots,from="symbol",
			   to="uniprot",species="mouse")
uniprot[1] <- "Q8VHW2"

# Loop to align proteins.
result <- list()
for (protB in uniprot) {
       result[[protB]] <- water(protA,protB)
}
names(result) <- paste(names(uniprot),uniprot,sep="|")

# Exctract scores.
scores <- sapply(result,function(x) x$Summary$Score)

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

