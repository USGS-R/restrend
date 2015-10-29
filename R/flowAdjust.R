#' Flow Adjusted Concentrations
#' 
#' Create flow adjusted concentrations for the seasonal Kendall 
#'trend test.
#' 
#' @param Stations a vector of the the station identifiers on which 
#'to do the flow adjustment.
#' @param Snames a vector of the response variables on which to do 
#'the flow adjustment.
#' @param use.logs fit a log-log LOWESS curve? If \code{FALSE}, then
#'do not use the log transforms.
#' @param span the span to use for the LOWESS curve.
#' @param max.cens the maximum percent censoring permitted. A warning
#'is printed if there are any censored data and a warning is printed
#'if any exceed the value.
#' @param report the name of the PDF file that contains graphs of 
#'each fit. The default is to use the name of the project with "_fa"
#'appended. If the PDF file exists, then it is not overwritten, 
#'but the name is appended with a sequence of numbers until one that
#'is valid is created.
#' @return The name of the report file.
#' @export
flowAdjust <- function(Stations="All", Snames="All", use.logs=TRUE,
											 span=0.75, max.cens=5, report) {
	## Coding history:
	##    2013Aug20 DLLorenz Original Coding
	##    2013Oct21 DLLorenz Finishing touches
	##
	## Preliminaries
	pos <- ckProj()
	Call <- match.call()
	if(missing(report)) {
		report <- paste(get("._Proj", pos=1), "_fa", sep="")
	}
	report.name <- setFileType(report, type="pdf")
	if(file.exists(report.name)) { # find one that works
		i <- 1L
		while(file.exists(report.name)) {
			report.tmp <- paste(report, zeropad(i, 2), sep="_")
			report.name <- setFileType(report.tmp, type="pdf")
			i <- i + 1L
		}
		report <- report.tmp
	}
	estrend.cl <- get("estrend.cl", pos=pos)
	estrend.in <- get("estrend.in", pos=pos)
	## Check to verify that data were set up for SK
	if(estrend.in$type == "tobit")
		stop("The project ", get("._Proj", pos=1), 
				 " was set up for the Tobit trend test.\n",
				 "No flow-adjustment necessary for the Tobit trend test.")
	estrend.df <- get("estrend.df", pos=pos)
	estrend.st <- get("estrend.st", pos=pos)
	estrend.cp <- get("estrend.cp", pos=pos)
	if(is.null(FLOW <- estrend.in$FLOW))
		stop("there is no FLOW defined for the current project")
	if(exists("estrend.fa", where=pos, inherits=FALSE)) {
		estrend.fa <- get("estrend.fa", pos=pos)
	} else {
		estrend.fa <- estrend.df # use as template
		for(station in estrend.in$stations)
			for(sname in estrend.in$snames)
				estrend.fa[station, sname][[1L]] <- c(0,0) # span and log
	}
	if(Stations[1L] == "All")
		Stations <- estrend.in$stations
	if(Snames[1L] == "All")
		Snames <- estrend.in$snames
	setPDF(basename=report)
	if(length(Snames) == 1L) {
		Ptest <- 1L
	} else if(length(Snames) == 2L) {
		par(mfrow=c(2,1))
		Ptest <- 2L
	} else if(length(Snames) <= 4L) {
		par(mfrow=c(2,2))
		Ptest <- 4L
	} else {
		par(mfrow=c(3,2))
		Ptest <- 6L
	} 
	on.exit(dev.off()) # close on error or normal exit
	## Track censoring
	Cens <- matrix("", nrow=0, ncol=2)
	Cens.pct <- double(0)
	## OK, do it by sname
	for(station in Stations) {
		i <- 0L
		for(sname in Snames) {
			if(estrend.st[station, sname] == "OK") {
				temp.df <- estrend.df[station, sname][[1L]]
				ckc <- 0
				if(censoring(temp.df[[sname]]) != "none") {
					Cens <- rbind(Cens, c(station, sname))
					Cens.pct <- c(Cens.pct, ckc <- estrend.cp[station, sname])
				}
				if(ckc <= max.cens) {
					Q <- temp.df[[FLOW]]
					Cn <- as.double(temp.df[[sname]])
					# Protect against 0s if log
					if(use.logs) {
						pick <- (Q > 0) & (Cn > 0)
						Q <- Q[pick]
						Cn <- Cn[pick]
					}
					AA.pl <- xyPlot(Q, Cn,
													Plot=list(size=0.05),
													xaxis.log=use.logs, yaxis.log=use.logs,
													xtitle="Streamflow", ytitle=sname,
													margin=c(4.1, 5.1, 2.5, 2.0))
					addTitle(station)
					AA.pl <- addSmooth(AA.pl, span=span, evaluation=101)
					## Now compute the FACs
					Fits <- interpLine(AA.pl, yfromx=temp.df[[FLOW]])
					if(use.logs) {
						FAC <- log(as.numeric(temp.df[[sname]])) - log(Fits)
					} else
						FAC <- as.numeric(temp.df[[sname]]) - Fits
					temp.df$FAC <- FAC
					estrend.df[station, sname][[1L]] <- temp.df
					estrend.fa[station, sname][[1L]] <- c(span, use.logs)
					i <- i + 1L
				}
			}
		}
		if((i <- i %% Ptest) != 0L)
			for(i in seq(Ptest - i + 1L))
				plot.new() # Fill up the page
	}
	## Finish up
	estrend.cl <- c(estrend.cl, Call)
	assign("estrend.df", estrend.df, 2)
	assign("estrend.fa", estrend.fa, 2)
	assign("estrend.cl", estrend.cl, 2)
	g.data.save()
	if(nrow(Cens))
		warning("Censored data in:\n",
						paste(paste(Cens[, 1L], Cens[, 2L], sep=", "), collapse="\n"))
	if(any(Cens.pct > max.cens))
		warning("Excessive censoring in:\n",
						paste(paste(Cens[Cens.pct > max.cens, 1L], 
												Cens[Cens.pct > max.cens, 2L], sep=", "), collapse="\n"))
	return(report.name)
}
