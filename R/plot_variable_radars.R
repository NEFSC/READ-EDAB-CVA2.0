#' @title Make Variable Importance Radars
#' @description
#' Produces radar plots of variable importance 
#'
#' @param variable_importance 
#' @param fig_name name to save figure as. Figure will save to current working directory if desired directory is not included in the figure name

#' @return Function does not return anything. Figures are saved to species-specific \code{figures} folder.
#'
#'@export
#'
plot_variable_radars <- function(variable_importance, fig_name){
  
  # Failsafe: If the function crashes, forcefully close any open PDFs
  on.exit(
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    },
    add = TRUE
  )
  
  variable_importance <- rbind(rep(1, ifelse(is.null(dim(variable_importance)), length(variable_importance), ncol(variable_importance))), 
                               rep(0, ifelse(is.null(dim(variable_importance)), length(variable_importance), ncol(variable_importance))), 
                               variable_importance)
                               
  if(nrow(variable_importance) == 3){
    col = 'grey'
  } else {
    col = 
  }                  
  
  grDevices::pdf(
    fig_name,
    width = 8,
    height = 8
  )     
  
  fmsb::radarchart(
    as.data.frame(variable_importance),
    pfcol = scales::alpha('grey', 0.5),
    seg = 10
  )
  
  grDevices::dev.off()
}