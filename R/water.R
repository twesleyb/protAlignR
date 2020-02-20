# Function to perform water alignment.
algnWater <- function(protA,protB){
	suppressPackageStartupMessages({
		library(seqinr)
	})
	# Download protein sequences in fasta format.
	baseurl <- "https://www.uniprot.org/uniprot/"
	urls <- paste0(baseurl,c(protA,protB),".fasta")
	invisible({
		sapply(urls,function(x) {
			       download.file(x,basename(x),quiet=TRUE) })
	})
	# Execute run-water to align proteins.
	fasta_files <- paste(basename(urls),collapse=" ")
	cmd <- paste("water",fasta_files,"-auto","-stdout")
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
