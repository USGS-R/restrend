#' Annual Trends
#' 
#' Perform the annual Kendall trend test.
#' 
#' The \code{kensen.test} includes a correction to the attained p-value
#'that accounts the first order auto regression of the data. It is 
#'applied whenever there are at least 10 observations. To suppress the
#'correction, set \code{ar.adj} to \code{FALSE}, or a numeric value can
#'be specified for the minimum number of observations, whihc can be 
#'useful for a ragged analysis. The value for \code{ar.adj} cannot be
#'set to less than 10.
#'
#' @param Stations a vector of the the station identifiers on which 
#'to do the flow adjustment.
#' @param Snames a vector of the response variables on which to do 
#'the flow adjustment.
#' @param use.logs log transform the data before the trend test?
#' @param ar.adj adjust the attained p-value to account for serial
#'correlation? See \bold{Details}.
#' @param report the name of the PDF file that contains a report for 
#'each test. The default is to use the name of the project with "_an"
#'appended. If the PDF file exists, then it is not overwritten, 
#'but the name is appended with a sequence of numbers until one that
#'is valid is created.
#' @return The name of the report file.
#' @export
annualTrends <- function(Stations="All", Snames="All", use.logs=FALSE,
											 ar.adj=TRUE, report) {
	## Coding history:
	##    2013Oct23 DLLorenz Original Coding
	##
	## Preliminaries
	pos <- ckProj()
	Call <- match.call()
	if(missing(report)) {
		report <- paste(get("._Proj", pos=1), "_an", sep="")
	}
	report.name <- setFileType(report, type="pdf")
	if(file.exists(report.name)) { # find one that works
		i <- 1L
		while(file.exists(report.name)) {
			report.tmp <- paste(report, zeroPad(i, 2), sep="_")
			report.name <- setFileType(report.tmp, type="pdf")
			i <- i + 1L
		}
		report <- report.tmp
	}
	if(ar.adj) {
		ar.adj <- max(10, ar.adj)
	} else
		ar.adj <- Inf
	estrend.cl <- get("estrend.cl", pos=pos)
	estrend.in <- get("estrend.in", pos=pos)
	## Check to verify that data were set up for annual
	if(estrend.in$type != "annual")
		stop("The project ", get("._Proj", pos=1), 
				 " was not set up for the annual Kendall trend test.")
	estrend.df <- get("estrend.df", pos=pos)
	estrend.st <- get("estrend.st", pos=pos)
	DoFAC <- FALSE
	if(!is.null(estrend.in$FLOW)) {
		if(exists("estrend.fa", where=pos)) {
			estrend.fa <- get("estrend.fa", pos=pos)
			DoFAC <- TRUE
		} else
			warning("Flow was defined for this project, but flow adjustment was not done")
	}
	DATES <- estrend.in$DATES
	DoReg <- estrend.in$analysis == "regular"
	## Get/create the output object
	if(exists("estrend.an", where=pos, inherits=FALSE)) {
		estrend.an <- get("estrend.an", pos=pos)
	} else {
		estrend.an <- estrend.df # use as template
		for(station in estrend.in$stations)
			for(sname in estrend.in$snames)
				estrend.an[station, sname][[1L]] <- list()
	}
	if(Stations[1L] == "All")
		Stations <- estrend.in$stations
	if(Snames[1L] == "All")
		Snames <- estrend.in$snames
	setPDF(basename=report)
	on.exit(dev.off()) # close on error or normal exit
	## Protect against warning errors in call to text
	warn <- options("warn")
	## OK, do it by sname
	for(station in Stations) {
		for(sname in Snames) {
			if(estrend.st[station, sname] == "OK") {
				temp.df <- estrend.df[station, sname][[1L]]
				## Convert to decimal time if necessary
				if(isDateLike(temp.df[[DATES]]))
					temp.df[[DATES]] <- dectime(temp.df[[DATES]])
				if(use.logs) {
					SKu <- eval(parse(text=paste("kensen.test(log(", sname, "),",
																			 DATES, ",", ar.adj, ")", sep="")),
											envir=temp.df)
				} else
					SKu <- eval(parse(text=paste("kensen.test(", sname, ",",
																			 DATES, ",", ar.adj, ")", sep="")),
											envir=temp.df)
				## Create report
				AA.lo <- setLayout(height=c(3,6))
				setGraph(1, AA.lo)
				par(mar=c(0,0,0,0), usr=c(0,1,0,1))
				txt <- c(paste(station, sname, sep="  "),
								 capture.output(print(SKu)))
				options(warn=-1)
				text(0, 1, paste(txt, collapse="\n"), family="mono", adj=c(0,1))
				options(warn)
				setGraph(2, AA.lo)
				AA.pl <- timePlot(temp.df[[DATES]], temp.df[[sname]], Plot=list(what="points"),
													ytitle=sname, yaxis.log=use.logs)
				refLine(coefficients=SKu$coef, current=AA.pl)
				ret <- list(aK=SKu)
				# OK, can we do FAC?
				if(DoFAC) {
					SKu <- eval(parse(text=paste("kensen.test(FAC,",
																			 DATES, ",", ar.adj,")", sep="")),
											envir=temp.df)
					## Create report
					AA.lo <- setLayout(height=c(3,6))
					setGraph(1, AA.lo)
					par(mar=c(0,0,0,0), usr=c(0,1,0,1))
					txt <- c(paste(station, sname, sep="  "),
									 capture.output(print(SKu)))
					options(warn=-1)
					text(0, 1, paste(txt, collapse="\n"), family="mono", adj=c(0,1))
					options(warn)
					setGraph(2, AA.lo)
					AA.pl <- timePlot(temp.df[[DATES]], temp.df[[sname]], Plot=list(what="points"),
														ytitle=sname, yaxis.log=use.logs)
					refLine(coefficients=SKu$coef, current=AA.pl)
					ret$aKFAC <- SKu
				} # End of FAC
				estrend.an[station, sname][[1L]] <- ret # Pack it up
			}
		}
	}
	## Finish up
	estrend.cl <- c(estrend.cl, Call)
	assign("estrend.cl", estrend.cl, pos=pos)
	assign("estrend.an", estrend.an, pos=pos)
	g.data.save()
	return(report.name)
}
