splitAt <- function(x, position) {
	return(unname(split(x, cumsum(seq_along(x) %in% position))))
}
