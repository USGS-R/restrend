#' Tobit Trends
#' 
#' Perform the tobit regression trend test.
#' 
#' @param Stations a vector of the station identifiers on which 
#'to do the flow adjustment.
#' @param Snames a vector of the response variables on which to do 
#'the flow adjustment.
#' @param use.logs log transform the data (and flow) before the trend test?
#' @param Flow logical, include the flow variable? Can also be an expression
#'that includes any transformation of flow.
#' @param Seasons include first-order fourier seasonal terms? Can also be
#'an integer value indicating the order of the terms to include.
#' @param Covars the covariates listed in the regression model.
#' @param report the name of the PDF file that contains a report for 
#'each test. The default is to use the name of the project with "_tb"
#'appended. If the PDF file exists, then it is not overwritten, 
#'but the name is appended with a sequence of numbers until one that
#'is valid is created.
#' @return The name of the report file.
#' @note The trend variable name is \code{.Dectime}, which represents
#'the trend in units of per year.
#' @export
tobitTrends <- function(Stations="All", Snames="All", use.logs=TRUE,
											 Flow=TRUE, Seasons=TRUE, Covars=NULL, report) {
	## Coding history:
	##    2013Sep17 DLLorenz Original Coding
	##    2013Oct21 DLLorenz Finishing touches
	##
	## Preliminaries
	FlowExp <- deparse(substitute(Flow))
	SeasonK <- as.integer(Seasons)
	pos <- ckProj()
	Call <- match.call()
	if(missing(report)) {
		report <- paste(get("._Proj", pos=1), "_tb", sep="")
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
	## Check to verify that data were set up for Tobit
	if(estrend.in$type != "tobit")
		stop("The project ", get("._Proj", pos=1), 
				 " was not set up for the Tobit trend test.")
	estrend.st <- get("estrend.st", pos=pos)
	estrend.df <- get("estrend.df", pos=pos)
	estrend.cn <- get("estrend.cn", pos=pos)
	estrend.cp <- get("estrend.cp", pos=pos)
	DATES <- estrend.in$DATES
	FLOW <- estrend.in$FLOW
	## Logic to express flow or log(flow) or any transform
	## FlowLog is logic to include flow, not log transform of flow
	FlowLog <- as.logical(FlowExp)
	if(is.na(FlowLog)) {
		FlowLog <- TRUE
	} else if(FlowLog) {
		if(use.logs) {
			FlowExp <- paste("log(", FLOW, ")", sep="")
		} else
			FlowExp <- FLOW
	}
	if(!is.null(Covars))
		Covars <- paste(Covars, collapse=" + ")
	## Get/create the output object
	if(exists("estrend.tb", where=pos, inherits=FALSE)) {
		estrend.tb <- get("estrend.tb", pos=pos)
	} else {
		estrend.tb <- estrend.df # use as template
		for(station in estrend.in$stations)
			for(sname in estrend.in$snames)
				estrend.tb[station, sname][[1L]] <- list()
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
			if(estrend.st[station, sname] == "OK" && 
				 	estrend.cp[station, sname] < 80) {
				temp.df <- estrend.df[station, sname][[1L]]
				temp.df$.Dectime <- dectime(temp.df[[DATES]])
				## Build the call
				if(use.logs) {
					form <- paste("log(", sname, ") ~ ", sep="")
				} else
					form <- paste(sname, " ~ ", sep="")
				if(FlowLog) {
						form <- paste(form, FlowExp, sep="")
				}
				if(Seasons) 
					form <- paste(form, " + fourier(.Dectime,", SeasonK, ")", sep="")
				if(!is.null(Covars))
					form <- paste(form, " + ", Covars, sep="")
				## Add decimal years and do it
				form <- paste("censReg(", form, " + .Dectime)", sep="")
				TB <- eval(parse(text=form), envir=temp.df)
				## Create report
				AA.lo <- setLayout(height=c(3,3.5,3.5), width=c(3.5,3.5))
				setGraph(1, AA.lo)
				par(mar=c(0,0,0,0), usr=c(0,1,0,1), xpd=NA)
				txt <- c(paste(station, sname, sep="  "),
								 capture.output(print(TB)))
				options(warn=-1)
				text(0, 1, paste(txt, collapse="\n"), family="mono", adj=c(0,1))
				options(warn)
				TBx <- summary(TB)
				par(xpd=FALSE)
				## The graphs
				## Overall fit
				setGraph(3, AA.lo)
				plot(TBx, which=2, set.up=FALSE)
				## Q-Normal plot for residuals
				setGraph(4, AA.lo)
				plot(TBx, which=5, set.up=FALSE)
				## The partial fit for time
				setGraph(5, AA.lo)
				plot(TBx, which=".Dectime", set.up=FALSE)
				## The partial fit for flow if possible
				setGraph(6, AA.lo)
				if(FlowLog) {
					plot(TBx, which=FlowExp, set.up=FALSE)
				} else
					plot.new()
				estrend.tb[station, sname][[1L]] <- list(TB=TB) # Pack it up
			}
		}
	}
	## Finish up
	estrend.cl <- c(estrend.cl, Call)
	assign("estrend.cl", estrend.cl, pos=pos)
	assign("estrend.tb", estrend.tb, pos=pos)
	g.data.save()
	return(report.name)
}
