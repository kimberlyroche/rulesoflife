#' Create output directory/ies that may not already exist
#'
#' @param path_to_dir a character vector of subdirectories to create (if not
#' already in existence)
#' @return the path as a joined string
#' @export
check_dir <- function(path_to_dir) {
  for(i in 1:length(path_to_dir)) {
    dir.create(do.call(file.path, as.list(path_to_dir[1:i])),
               showWarnings = FALSE)
  }
  return(do.call(file.path, as.list(path_to_dir)))
}
