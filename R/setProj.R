#' Set Up an ESTREND Project
#' 
#' Define the data and other characteristics of an ESTREND project.
#' 
#' If \code{Start} and \code{End} are \code{NULL}, then a ragged 
#'analysis is set up, and each variable set up for each station is 
#'analyzed on the available data. Otherwise, the analysis is regular and
#'will be restricted to the time frame specified by \code{Start} and 
#'\code{End}.
#'
#' If \code{type} is "seasonal," then the data are processed for a 
#'seasonal Kendall type of analysis---seasons are defined (12, 6, 4, 
#'and 3 per year) and evaluated to select the "best" number of seasons.
#'If \code{type} is "tobit" or "annual," then no seasonal analysis is
#'done because the data are ready for analysis. If \code{type} is
#'"annual," then the data must be uncensored.
#'
#' @param project the name of the project to set up. The project name 
#'is forced to all lower case. 
#' @param data the dataset to use in for the project.
#' @param STAID the name of the column in \code{data} that contains the 
#'station identifiers. The data are forced to character for the analysis.
#' @param DATES the name of the column in \code{data} containing the 
#'sample dates. This column must be class "Date" for seasonal analyses,
#'but can be numeric or integer for annual analyses.
#' @param Snames the name of the columns in \code{data} containing the 
#'sample data for trend analysis. These must be of class "numeric,"
#'"integer," "lcens," or "qw."
#' @param FLOW the name of the column in \code{data} containing the 
#'streamflow for each sample.
#' @param Covars the name of the columns in \code{data} containing any 
#'covariate data for trend analysis.
#' @param type the kind of analysis. Must be "seasonal," "tobit," or
#'"annual." Only the fist letter is necessary. See \bold{Details}.
#' @param Start the starting date for the analysis. For seasonal analyses,
#'must be "Date" or a character string the represents a date. For
#'annual analyses, must match the type of \code{DATES}. 
#'See \bold{Details}.
#' @param End the ending date for the analysis. For seasonal analyses,
#'must be "Date" or a character string the represents a date. For
#'annual analyses, must match the type of \code{DATES}. To guarantee
#'that the periods are set up correctly, \code{End} should be the first
#'day of the month following the actual last day.
#' @param tol the tolerance for the samples dates. To be included in a
#'regular analysis, the first sample within the \code{Start} to 
#'\code{End} time period must be within \code{tol} and likewise for the
#'last sample. if \code{tol} is \code{NULL}, then it is set to 5 percent
#'of the time frame, except for annual series analysis when it is set to
#'1 year.
#' @param min.obs the minimum number of observations required for a 
#'trend analysis.
#' @return The name of the project is returned.
#' @note A directory is created using the name of the project is created
#'in the user's directory. It contains the objects created by restrend
#'as R workspaces.
#' @export
setProj <- function(project, data, STAID, DATES, Snames, FLOW=NULL,
										Covars=NULL, type="seasonal", Start=NULL,
										End=NULL, tol=NULL, min.obs=20) {
	## Coding history:
	##    2013Aug20 DLLorenz Original Coding
	##    2013Oct21 DLLorenz Finishing touches
	##
	## Test the project
	project <- tolower(project)
	if(file.exists(project))
		stop("A directory or file named ", project, 
				 " exists in the user's directory.")
	## Check search path for a current project to detach it
	if(exists("._Proj", where=1L)) {
		cur <- get("._Proj", pos=1L)
		if(cur %in% search())
			detach(cur, character.only=TRUE)
	}
	## Create the project
	g.data.attach(project, warn=FALSE)
	on.exit(detach(2)) # If error, trap false creation
	## Trim dataset if necessary
	analysis <- "ragged"
	if(!is.null(c(Start, End))) { # Both not null
		if(is.null(Start) || is.null(End)) 
			stop("Either both or neither Start and End must be specified.")
		## Fix to allow End to be the end of the data for annual analysis.
		if(type == "annual" && !isDateLike(data[[DATES]])) {
			data <- data[data[[DATES]] >= Start & data[[DATES]] <= End, ]
		} else {
			## Check for consistent begin/end when "seasonal"
			if(type == "seasonal") {
				if(month(as.POSIXlt(Start)) != month(as.POSIXlt(End)) ||
					 	day(as.POSIXlt(Start)) != day(as.POSIXlt(End)))
					stop("The Start and End day and month must agree for seasonal analyses")
			}
			data <- data[data[[DATES]] >= Start & data[[DATES]] < End, ]
		}
		analysis <- "regular"
		if(is.null(tol)) {
			if(is.character(Start))
				Start <- as.Date(Start)
			if(is.character(End))
				End <- as.Date(End)
			tol <- (End - Start)/20
			if(type == "annual" && !isDateLike(Start))
				tol <- 1
		}
	}
	## Create the data structures
	Call <- match.call()
	if(any(bad <- !(Snames %in% names(data))))
		stop(paste(Snames[bad], collapse=", "), " not found in data")
	data[[STAID]] <- as.character(data[[STAID]])
	stations <- unique(data[[STAID]])
	estrend.by <- expand.grid(stations, Snames)
	names(estrend.by) <- c("stations", "snames")
	estrend.df <- by(estrend.by, estrend.by, as.data.frame)
	for(station in stations)
		for(col in Snames) {
			estrend.df[station, col][[1L]] <- 
				na.omit(data[data[[STAID]] == station, 
										 c(STAID, DATES, col, FLOW, Covars)])
		}
	## Now the status and censoring info
	estrend.st <- matrix("", nrow=length(stations), ncol=length(Snames),
											 dimnames=list(stations=stations, snames=Snames))
	estrend.cn <- estrend.st # Just copy
	estrend.cp <- matrix(0, nrow=length(stations), ncol=length(Snames),
											 dimnames=list(stations=stations, snames=Snames))
	for(station in stations)
		for(col in Snames) {
			temp.df <- estrend.df[station, col][[1]]
			if(nrow(temp.df) == 0) {
				estrend.st[station, col] <- "no data"
				estrend.cn[station, col] <- "none"
			} else if(nrow(temp.df) < min.obs) {
				estrend.st[station, col] <- "too few data"
				estrend.cn[station, col] <- censoring(temp.df[[col]])
			} else if(analysis == "regular" && max(temp.df[[DATES]]) < End - tol) {
				estrend.st[station, col] <- "short record"
				estrend.cn[station, col] <- censoring(temp.df[[col]])
			} else if(analysis == "regular" && min(temp.df[[DATES]]) > Start + tol) {
				estrend.st[station, col] <- "short record"
				estrend.cn[station, col] <- censoring(temp.df[[col]])
			} else {
				estrend.st[station, col] <- "OK"
				estrend.cn[station, col] <- censoring(temp.df[[col]])
			}
			## Compute the total percentage censoring
			if(estrend.cn[station, col] != "none")
				estrend.cp[station, col] <- pctCens(temp.df[[col]])
		}
	## And the analysis info
	type=match.arg(type, c("seasonal", "tobit", "annual"))
	estrend.in <- list(stations=stations, snames=Snames,
										 analysis=analysis, type=type,
										 Start=Start, End=End, DATES=DATES,
										 FLOW=FLOW, Covars=Covars)
	## Save the data
	assign("estrend.df", estrend.df, 2)
	assign("estrend.st", estrend.st, 2)
	assign("estrend.cn", estrend.cn, 2)
	assign("estrend.cp", estrend.cp, 2)
	assign("estrend.in", estrend.in, 2)
	## Process the seasons if needed.
	if(type == "seasonal") {
		## Season listing and season select
		estrend.sl <- by(estrend.by, estrend.by, function(x) list())
		estrend.ss <- matrix(0, nrow=length(stations), ncol=length(Snames),
												 dimnames=list(stations=stations, snames=Snames))
		for(station in stations)
			for(col in Snames) {
				temp.df <- estrend.df[station, col][[1]]
				if(estrend.st[station, col] == "OK") {
					if(analysis == "regular") {
						ssn12 <- regularSeries(seq(nrow(temp.df)), 
																	 temp.df[[DATES]], 
																	 begin=Start,
																	 end=End, k.period=1)$Value
					} else {
						ssn12 <- regularSeries(seq(nrow(temp.df)), 
																	 temp.df[[DATES]], k.period=1)$Value
						if(ck12 <- (length(ssn12) %/% 12) > 0)
							ssn12 <- c(ssn12, rep(NA_real_, 12L - ck12))
					}
					temp.l <- estrendSeasonTable(ssn12)
					estrend.sl[station, col][[1]] <- temp.l
					estrend.ss[station, col][[1]] <- temp.l$best
				}
			}
		## Save the results
		assign("estrend.sl", estrend.sl, 2)
		assign("estrend.ss", estrend.ss, 2)
	}
	# Preserve a history of this project.
	assign("estrend.cl", list(Call), 2) 
	g.data.save()
	assign("._Proj", project, pos=1L)
	## Recover and exit
	on.exit()
	return(project)
}
