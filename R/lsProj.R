#' List ESTREND Projects
#' 
#' List the ESTREND projects that have been set up in the current workspace.
#' 
#' @return The names of valid projects are returned. If a project is in use,
#'then it is labeled "in use," otherwise the most recently used project is
#'labeled "recent" if that can be determined.
#' @export
lsProj <- function() {
	## Coding history:
	##    2013Oct22 DLLorenz Original Coding
	##
  dirs <- list.dirs(".", full.names=FALSE, recursive=FALSE)
  ## Strip leading "./"
  dirs <- sub("./", "", dirs, fixed=TRUE)
  projs <- paste(dirs, "/estrend.in.RData", sep="")
  projs <- file.exists(projs)
	## Any in use or most recent?
  retval <- dirs[projs]
  if(length(retval)) {
  	names(retval) <- rep("", length(retval))
  	if(exists("._Proj", where=1L)) {
  		cur <- get("._Proj")
  		if(cur %in% search()) {
  			names(retval)[retval == cur] <- "in use"
  		} else
  			names(retval)[retval == cur] <- "recent"
  	}
  } else
  	retval <- c(Error="No projects have been set")
	return(retval)
}
