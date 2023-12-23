
convert_gt <- function(gt) {
  if (gt == "0/0") {
    return(0)
  } else if (gt == "0/1" || gt == "1/0") {
    return(1)
  } else if (gt == "1/1") {
    return(2)
  } else {
    return(NA)  # NA for missing or non-diploid genotypes
  }
}

