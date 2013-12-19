#' Plot Trends
#' 
#' Create a series of diagnostic plots for a single Tobit trend test.
#' 
#' @param data the data set returned from \code{getTrends}.
#' @param which which trend to plot, either "Trend" or "Trend.pct."
#' @param device the name of the graphics device. Defaults to a valid
#'graphics device on any platform. May be "pdf" to create a pdf file.
#' @return The name of the graphics device.
#' @note The project from which \code{data} was retrieved must be the 
#'current project.
#' @export
plotTrends <- function(data, which="Trend", device) {
	## Coding history:
	##    2013Oct21 DLLorenz Original Coding
	##
	## Preliminaries
	Stations <- unique(data$Station)
	Snames <- unique(data$Response)
	which <- match.arg(which, c("Trend", "Trend.pct"))
	device.name <- "Trends"
	pos <- ckProj()
	estrend.df <- get("estrend.df", pos=pos)
	estrend.in <- get("estrend.in", pos=pos)
	DATES <- estrend.in$DATES
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
	for(station in Stations) {
		for(sname in Snames) {
			temp.df <- estrend.df[station, sname][[1L]]
			## Plot the data
			AA.pl <- timePlot(temp.df[[DATES]], temp.df[[sname]],
												Plot=list(name="Observed data", what="points", size=0.07),
												ytitle=sname)
			## Plot the trend(s)
			Dates <- temp.df[[DATES]]
			DecDates <- dectime(Dates)
			midDate <- mean(range(Dates))
			RelDates <- DecDates - dectime(midDate)
			Tnd <- data[data$Station == station & data$Response == sname, which]
			Type <- data[data$Station == station & data$Response == sname, "Type"]
			RepValue <- data[data$Station == station & data$Response == sname, "RepValue"]
			Colors <- blueRed.colors(length(Type)  + 1L)
			Dirs <- data[data$Station == station & data$Response == sname, "Trend.dir"]
			for(i in seq(length(Type))) {
				if(which == "Trend") {
					y <- Tnd[i]*RelDates + RepValue[i]
				} else
					y <- exp(log(Tnd[i]/100 + 1)*RelDates)*RepValue[i]
				## set bold line if significant
				if(Dirs[i] %in% c("up", "down")) {
					Wid <- "bold"
				} else
					Wid <- "standard"
				AA.pl <- addXY(Dates, y, Plot=list(name=Type[i], what="lines", 
																					 color=Colors[i], width=Wid),
											 current=AA.pl)
			}
			## Place the explanation
			if(Tnd[1L] > 0) {
				addExplanation(AA.pl, where="ul", title="")
			} else
				addExplanation(AA.pl, where="ur", title="")
		}
	}
	if(closeIt)
		dev.off()
	par(askIt)
	return(device.name)
}
