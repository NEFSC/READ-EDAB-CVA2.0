#' @title Helper: Force Data.Frame to Wide
#' @description
#' Forces Input Data into Wide format for summary tables
#'
#' @param df data.frame containing data to force into wide format
#' @param raw_nms a character string of raw variable names to match to
#' @param sc name of column indicating stocks to match to data
#' @param stocks_obj shpfile of sttocks. used to match stock names to data
#'
#' @return a wide data.frame
#'
#'@export
#'

# Forces data into wide format no matter how it arrives (vertical, long, or named vector)
force_wide <- function(df, raw_nms, sc, stocks_obj) {
  if (is.numeric(df) && !is.matrix(df)) {
    df <- as.data.frame(t(df), check.names = FALSE)
  }
  df <- as.data.frame(df, check.names = FALSE)

  colnames(df) <- repair_names(colnames(df), raw_nms)
  rownames(df) <- repair_names(rownames(df), raw_nms)

  # Case 1: Already Wide
  if (any(raw_nms %in% colnames(df))) {
    # Do nothing
    # Case 2: Vertical Format (terra::global drops variables into row names)
  } else if (any(raw_nms %in% rownames(df))) {
    df <- as.data.frame(t(df), check.names = FALSE)

    # Case 3: Long Format (variables are in a column)
  } else {
    name_col <- NULL
    for (c in colnames(df)) {
      if (is.character(df[[c]]) || is.factor(df[[c]])) {
        if (any(raw_nms %in% as.character(df[[c]]))) {
          name_col <- c
          break
        }
      }
    }
    if (!is.null(name_col)) {
      num_col <- names(which(sapply(df, is.numeric)))[1]
      if (!is.null(num_col)) {
        df <- tidyr::pivot_wider(
          df,
          names_from = dplyr::all_of(name_col),
          values_from = dplyr::all_of(num_col)
        )
      }
    }
  }

  # Inject missing stock columns
  if (!sc %in% colnames(df)) {
    df[[sc]] <- if (!is.null(stocks_obj) && !is.null(stocks_obj$stock_area)) {
      stocks_obj$stock_area[1]
    } else {
      "global"
    }
  }
  return(df)
}
