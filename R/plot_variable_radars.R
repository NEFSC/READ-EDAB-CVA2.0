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
  
  if(is.null(dim(variable_importance))){
    np = 1
  } else {
    np = nrow(variable_importance)
  }
  
  variable_importance <- rbind(rep(1, ifelse(is.null(dim(variable_importance)), length(variable_importance), ncol(variable_importance))), 
                               rep(0, ifelse(is.null(dim(variable_importance)), length(variable_importance), ncol(variable_importance))), 
                               variable_importance)
                               
  if(np == 1){
    col = 'grey'
  } else {
    col = RColorBrewer::brewer.pal(n = np, 'Set1')
  }                  
  
  grDevices::pdf(
    fig_name,
    width = 8,
    height = 8
  )     
  
  fmsb::radarchart(
    as.data.frame(variable_importance),
    pfcol = scales::alpha(col,0.2),
    seg = 10,
    pty = 15:(15+(np-1)),
    pcol = scales::alpha(col, 1),
    plty = 1,
    vlcex = 1.25
  )
  
  if(np > 1){
    graphics::legend(
      'bottomright',
      legend = rownames(variable_importance)[3:nrow(variable_importance)],
      lty = 1,
      col = col,
      pch = 15:(15+(np-1)),
      bty = 'n'
    )
  }
  
  grDevices::dev.off()
}