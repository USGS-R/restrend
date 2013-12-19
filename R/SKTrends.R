#' Seasonal Kendall Trends
#' 
#' Perform the seasonal Kendall trend test.
#' 
#' @param Stations a vector of the the station identifiers on which 
#'to do the flow adjustment.
#' @param Snames a vector of the response variables on which to do 
#'the flow adjustment.
#' @param use.logs log transform the data before the trend test, 
#'applies only to uncensored seasonal Kendall test?
#' @param max.cens the maximum percent censoring permitted for the
#'uncensored seasonal Kendall test. If the percentage od censoring
#'exceeds this value, then the censored seasonal Kendall test is
#'performed. Set to a negative value to force the censored seasonal 
#'Kendall test for all \code{Stations} and \code{Snames}.
#' @param nseas the number of seasons to use for all of the tests. If
#'\code{NULL} (default), then use the selected number of seasons 
#'defined in \code{setProj}.
#' @param report the name of the PDF file that contains a report for 
#'each test. The default is to use the name of the project with "_sk"
#'appended. If the PDF file exists, then it is not overwritten, 
#'but the name is appended with a sequence of numbers until one that
#'is valid is created.
#' @return The name of the report file.
#' @export
SKTrends <- function(Stations="All", Snames="All", use.logs=TRUE,
											 max.cens=5, nseas=NULL, report) {
	## Coding history:
	##    2013Sep05 DLLorenz Original Coding
	##
	## Preliminaries
	pos <- ckProj()
	Call <- match.call()
	if(missing(report)) {
		report <- paste(get("._Proj", pos=1), "_sk", sep="")
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
	if(estrend.in$type != "seasonal")
		stop("The project ", get("._Proj", pos=1), 
				 " was not set up for the seasonal Kendall trend test.")
	estrend.df <- get("estrend.df", pos=pos)
	estrend.st <- get("estrend.st", pos=pos)
	estrend.ss <- get("estrend.ss", pos=pos)
	estrend.cn <- get("estrend.cn", pos=pos)
	estrend.cp <- get("estrend.cp", pos=pos)
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
	if(exists("estrend.sk", where=pos, inherits=FALSE)) {
		estrend.sk <- get("estrend.sk", pos=pos)
	} else {
		estrend.sk <- estrend.df # use as template
		for(station in estrend.in$stations)
			for(sname in estrend.in$snames)
				estrend.sk[station, sname][[1L]] <- list()
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
				if(is.null(nseas)) {
					NumSeas <- estrend.ss[station, sname]
				} else
					NumSeas <- nseas
				ckc <- estrend.cp[station, sname]
				if(ckc <= max.cens) {
					if(DoReg) {
						RegSer <- regularSeries(as.double(temp.df[[sname]]), 
																		temp.df[[DATES]],
																		begin=estrend.in$Start,
																		end=estrend.in$End, 
																		k.period=12/NumSeas)
					} else
						RegSer <- regularSeries(as.double(temp.df[[sname]]), 
																		temp.df[[DATES]],
																		k.period=12/NumSeas)
					## Do it, finagles used to maintain names
					names(RegSer)[4L] <- sname
					if(use.logs) {
						SKu <- eval(parse(text=paste("seaken(log(", sname, "),",
																						NumSeas, ")", sep="")),
												envir=RegSer)
					} else
						SKu <- eval(parse(text=paste("seaken(", sname, ",",
																						NumSeas, ")", sep="")),
												envir=RegSer)
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
					seriesPlot(SKu)
					ret <- list(SKu=SKu)
					# OK, can we do FAC?
					if(DoFAC) {
						if(DoReg) {
							RegSer <- regularSeries(as.double(temp.df[["FAC"]]), 
																			temp.df[[DATES]],
																			begin=estrend.in$Start,
																			end=estrend.in$End, 
																			k.period=12/NumSeas)
						} else
							RegSer <- regularSeries(as.double(temp.df[["FAC"]]), 
																			temp.df[[DATES]],
																			k.period=12/NumSeas)
						## Do it, maintain names
						Name <- paste("FAC", sname, sep=".")
						names(RegSer)[4L] <- Name
						SKu <- eval(parse(text=paste("seaken(", Name, ",",
																						NumSeas, ")", sep="")),
												envir=RegSer)
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
						seriesPlot(SKu)
						ret$SKFAC <- SKu
					} # End of FAC
				} else { # do not meet censoring criterion
					temp.df$.ndx. <- seq(nrow(temp.df)) # Index to censored data
					if(DoReg) {
						RegSer <- regularSeries(temp.df$.ndx., 
																		temp.df[[DATES]],
																		begin=estrend.in$Start,
																		end=estrend.in$End, 
																		k.period=12/NumSeas)
					} else
						RegSer <- regularSeries(temp.df$.ndx., 
																		temp.df[[DATES]],
																		k.period=12/NumSeas)
					## Replace Value with the censored data
					RegSer$Value <- as.lcens(temp.df[[sname]])[RegSer$Value]
					## Do it, finagles used to maintain names
					names(RegSer)[4L] <- sname
					SKu <- eval(parse(text=paste("censSeaken(", sname, ",",
																			 NumSeas, ")", sep="")),
											envir=RegSer)
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
					seriesPlot(SKu)
					ret <- list(SKc=SKu)
				}
				estrend.sk[station, sname][[1L]] <- ret # Pack it up
			}
		}
	}
	## Finish up
	estrend.cl <- c(estrend.cl, Call)
	assign("estrend.cl", estrend.cl, pos=pos)
	assign("estrend.sk", estrend.sk, pos=pos)
	g.data.save()
	return(report.name)
}
