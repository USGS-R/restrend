#' Summarize Samples
#' 
#' Create a multi-page pdf file of sample data by station in a dataset.
#'The first pages are a listing of the first and last sample and total
#'numerof of samples at each station. The following pages are dot plots 
#'of the sample dates. No more than 40 stations per page are liested or 
#'plotted.
#' 
#' @param data the dataset to summarize
#' @param DATES the name of the column containing the sample dates
#' @param STAID the name of the column containing the station identifiers
#' @param file the output file base name; the .pdf suffix 
#'is appended to make the actual file name. If missing, then the
#'name of \code{data} is used as the base name.
#' @return The actual file name is returned invisibly.
#' @export
sampReport <- function(data, DATES="DATES", STAID="STAID", file) {
	## Coding history:
	##    2013Aug19 DLLorenz Original Coding
	##
	if(missing(file))
		file <- deparse(substitute(data))
	retval <- setPDF(basename=file)
	## Compute the summary stats
	summ <- sumSamp(data, DATES=DATES, STAID=STAID, by.numeric=FALSE)
	## Set up and Draw the text
	N <- nrow(summ)
	Pages <- (N - 1) %/% 40 + 1
	Ngrp <- N %/% Pages + 1
	Grps <- seq(N) %/% Ngrp
	for(i in unique(Grps)) {
		plot.new()
		par(mar=c(0,0,0,0), usr=c(0,1,0,1))
		txt <- capture.output(print(summ[Grps == i, ]))
		text(0, 1, paste(txt, collapse="\n"), family="mono", adj=c(0,1))
	}
	## Now the dot plots
	for(i in unique(Grps)) {
		subdata <- data[data[[STAID]] %in% (summ$STAID[Grps == i]),]
		dotPlot(subdata[[DATES]], subdata[[STAID]])
	}
	## All done, close the graph
	dev.off(retval[[1]])
	invisible(paste(retval[[2]], ".pdf", sep=""))
}
