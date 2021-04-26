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

#' Return a human-readable taxon label
#'
#' @param tax taxonomy data.frame
#' @param coord index of logratio to label
#' @param coord_system_label one of c("CLR", "ALR")
#' @return label
#' @export
get_tax_label <- function(tax, coord, coord_system_label) {
  coord_system_label <- toupper(coord_system_label)
  if(coord_system_label %in% c("CLR", "ALR")) {
    if(coord > nrow(tax)) {
      stop("Invalid taxonomy index")
    }
    tax_pieces <- tax[coord,]
    max_level <- max(which(!is.na(tax_pieces)))
    if(coord_system_label == "CLR") {
      label <- paste0(coord_system_label, "(ASVs in ", names(tax)[max_level], " ", tax_pieces[max_level], ")")
    } else {
      max_level_ref <- max(which(!is.na(tax[nrow(tax),])))
      # ALR
      numerator <- paste0("ASVs in ", names(tax)[max_level], " ", tax_pieces[max_level])
      denominator <- paste0("ASVs in ", names(tax)[max_level_ref], " ", tax[nrow(tax),max_level])
      label <- paste0("ALR(", numerator, "/", denominator, ")")
    }
    return(label)
  } else {
    stop("Unrecognized coordinate system!")
  }
}

#' Calculate a score (0.0 to 1.0) for the "universality" of given association
#' pair
#'
#' @param x vector of correlations between a pair of CLR microbes across hosts
#' @return numeric "score"
#' @export
calc_universality_score <- function(x) {
  x.sign <- sapply(x, sign)
  neg.idx <- which(x.sign < 0)
  pos.idx <- which(x.sign > 0)
  neg.abs.median <- NA
  pos.abs.median <- NA
  if(length(neg.idx) > 0) {
    neg.abs.median <- abs(median(x[neg.idx]))
  }
  if(length(pos.idx) > 0) {
    pos.abs.median <- abs(median(x[pos.idx]))
  }
  score <- 0
  # R doesn't do short circuit evaluation? yes, it does: update with &&
  if(!is.na(neg.abs.median) & !is.na(pos.abs.median)) {
    if(length(pos.idx) > length(neg.idx)) {
      # use this as majority direction
      score <- (length(pos.idx)/length(x))*(pos.abs.median)
    } else {
      score <- (length(neg.idx)/length(x))*(neg.abs.median)
    }
  } else {
    if(is.na(neg.abs.median)) {
      score <- (length(pos.idx)/length(x))*(pos.abs.median)
    } else {
      score <- (length(neg.idx)/length(x))*(neg.abs.median)
    }
  }
  score
}
