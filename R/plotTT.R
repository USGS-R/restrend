#' Plot Tobit Trends
#' 
#' Create a series of diagnostic plots for a single Tobit trend test.
#' 
#' @param Station the station identifier.
#' @param Sname the response variable.
#' @param device the name of the graphics device. Defaults to a valid
#'graphics device on any platform. May be "pdf" to create a pdf file.
#' @return The name of the graphics device.
#' @seealso \code{\link{plot.summary.censReg}}
#' @export
plotTT <- function(Station, Sname, device) {
	## Coding history:
	##    2013Oct18 DLLorenz Original Coding
	##
	## Preliminaries
	Station
	Sname
	device.name <- make.names(paste(Station, Sname, sep="_"))
	pos <- ckProj()
	estrend.tb <- get("estrend.tb", pos=pos)
	closeIt <- FALSE
	askIt <- par("ask")
	if(missing(device)) {
		setGD(device.name)
	} else if(device == "pdf") {
		setPDF(basename=device.name)
		device.name <- setFileType(device.name, type="pdf")
		closeIt <- TRUE
	} else {
		setPage(layout=list(width=6, height=6),name=device.name, device=device)
		par(ask=TRUE)
	}
	plot(summary(estrend.tb[Station, Sname][[1]]$TB), set.up=FALSE)
	if(closeIt)
		dev.off()
	par(askIt)
	return(device.name)
}
