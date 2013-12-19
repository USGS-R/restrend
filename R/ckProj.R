#' Check the ESTREND Project
#' 
#' Check and get the ESTREND project to use for this analyses.
#' 
#'
#' @return The position of the current project in the serach path 
#'is returned.
#' @export
ckProj <- function() {
	## Coding history:
	##    2013Aug20 DLLorenz Original Coding
	##
	## Test/Open the project
	if(exists("._Proj", where=1)) {
		project <- get("._Proj", pos=1)
	} else
		stop("No defined project.")
	pos <- which(search() == project)
	if(length(pos) == 0L)
		stop("No current project.")
	return(pos)
}
