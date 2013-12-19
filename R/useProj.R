#' Define the ESTREND Project
#' 
#' Define the ESTREND project to use for subsequent analyses.
#' 
#'
#' @param project the name of the project to use. The project name 
#'is forced to all lower case. If missing, then the most recent project
#'is opened for use.
#' @return The name of the project is returned.
#' @export
useProj <- function(project) {
	## Coding history:
	##    2013Aug20 DLLorenz Original Coding
	##    2013Oct21 DLLorenz Try to close existing project
	##
	## Test/Open the project
	if(missing(project))
		project <- get("._Proj", pos=1L)
	project <- tolower(project)
	if(!file.exists(project))
		stop("A directory or file named ", project, 
				 " does not exist in the user's directory.")
	## Check search path
	if(exists("._Proj", where=1L)) {
		cur <- get("._Proj", pos=1L)
		if(cur %in% search())
			detach(cur, character.only=TRUE)
	}
	if(project %in% search())
		detach(project, character.only=TRUE)
	g.data.attach(project) # Always attached to position 2
	if(!exists("estrend.cl", where=2L)) {
		detach(2L)
		stop(project, " is not a valid project.")
	}
	assign("._Proj", project, pos=1L)
	## Exit
	return(project)
}
