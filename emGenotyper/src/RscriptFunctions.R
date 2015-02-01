check.package.installation <- function(package.name, cran.repos, version=0, 
                                       dependencies=TRUE, quietly=TRUE){
  if(missing(package.name)) stop("No package name was supplied")
  if(missing(cran.repos)) {
    warning(paste("No CRAN repository mirror was specified, setting from", 
                  "getOption(\"repos\").  If none is set, you may be",
                  "prompted to choose a CRAN mirror.", sep=" "))
    cran.repos <- getOption("repos")
  }
  if(cran.repos == "@CRAN@") stop("Failed to set a CRAN mirror")
  
  if( ! is.element(package.name, installed.packages()[, 1]) ) {
    # if package doesn't exist at all, install it
    cat("Cannot find ", package.name, ", installing now\n", sep="")
    install.packages(package.name, repos=cran.repos, quiet=quietly,
                     dependencies=dependencies)
  } else if(packageVersion(package.name) < version) {
    # if package is outdated, installs the latest version
    cat("Can't find version >= ", version, " of ", package.name, 
        ", installing now\n", sep="")
    install.packages(package.name, repos=cran.repos, quiet=quietly,
                     dependencies=dependencies)
  }
  
  # load the package
  cat("Loading ", package.name, "\n", sep="")
  require(package.name,character.only=TRUE,quietly=quietly)
}
