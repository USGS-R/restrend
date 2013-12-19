#' Flow Adjusted Concentrations
#' 
#' Remove flow adjusted concentrations for the seasonal Kendall 
#'trend test.
#' 
#' The arguments for \code{Stations} and \code{Snames} must be valid
#'station identifiers and response varaibles. "All" is not valid.
#' @param Stations a vector of the the station identifiers on which 
#'to do the flow adjustment.
#' @param Snames a vector of the response variables on which to do 
#'the flow adjustment.
#' @return Nothing is returned.
#' @note The \code{undoFA} function is used to remove the flow adjusted
#'concentration data and prevent that analysis for the uncensored seasonal
#'Kendall test. It is intended for the rare case when a reasonable 
#'flow-adjustment model cannot be found.
#' @export
undoFA <- function(Stations, Snames) {
	## Coding history:
	##    2013Aug21 DLLorenz Original Coding
	##    2013Oct21 DLLorenz Finishing touches
	##
	## Preliminaries
	pos <- ckProj()
	Call <- match.call()
	estrend.cl <- get("estrend.cl", pos=pos)
	estrend.fa <- get("estrend.fa", pos=pos)
	estrend.df <- get("estrend.df", pos=pos)
	## OK, do it by sname
	for(station in Stations) {
		for(sname in Snames) {
			temp.df <- estrend.df[station, sname][[1L]]
			temp.df$FAC <- NULL
			estrend.df[station, sname][[1L]] <- temp.df
			estrend.fa[station, sname][[1L]] <- c(0, 0)
		}
	}
	## Finish up
	estrend.cl <- c(estrend.cl, Call)
	assign("estrend.df", estrend.df, 2)
	assign("estrend.fa", estrend.fa, 2)
	assign("estrend.cl", estrend.cl, 2)
	g.data.save()
	invisible()
}
