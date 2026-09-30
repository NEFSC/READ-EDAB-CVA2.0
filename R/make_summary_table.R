#' @title Make Vulnerability Summary Table
#' @description
#' Generates the summary vulnerability table with final factor scores, standard deviations, and distribution of spatial scores for each species, as well as final exposure scores. One table is generated per species and saved in the species folder. Requires directory to be set up per directions in the package documentation/manual.
#'
#' @param species names of the species to plot. Must match folder name to pull correct data.
#' @param forecast_release,hindcast_release MOM6 release codes for the (f)orecast and (h)indcasts used. Used to pull correct variable exposures
#' @param forecast_init forecast_initialization code corresponding to the forecast_initalization date of the desired forecast data. Used to pull correct variable exposures
#' @param hindcast_yr_range character string corresponding to the years in the hindcast data used. Used to pull correct ranked exposure values and save the data properly
#' @param stock_key a named vector containing the abbreviations and long names of stocks to use as labels
#' @param stock_order a vector containing the long stock names in the desired order
#' @param table_dir file path to folder to save tables in
#'
#' @return Function does not return anything. The table is saved as a png to the \code{table_dir} folder.
#'
#'@export
#'
make_summary_table <- function(species,
    metric,
    mean_metric,
    certainty_metric,
    raw_data,
    raw_names,
    clean_names,
    stocks,
    stock_key,
    stock_order,
    stock_col,
    certainty_col,
    total_sens_col,
    var_imp = NULL,
    imp_threshold = 0.1,
    save = T, 
    table_dir
    ) {
  
  attr_map <- setNames(clean_names, raw_names)
  
  #extract raw values from raster if metric is exposure or vulnerability
  if(metric != 'sensitivity'){
    # 1. Extract raster pixels within each stock polygon
    if (!is.null(stocks)) {
      # This returns a dataframe with 'ID' (polygon index) and one column per raster layer
      ext_vals <- terra::extract(raw_data, stocks)
      
      # Map the polygon ID back to the actual stock names
      ext_vals$stock <- stocks$stock_area[ext_vals$ID]
      ext_vals$ID <- NULL
      
      #extract all pixels
      global_vals <- as.data.frame(terra::values(raw_data))
      global_vals$stock <- 'global'
      
      all_vals <- dplyr::bind_rows(ext_vals, global_vals)
    } else {
      all_vals <- as.data.frame(terra::values(raw_data))
      all_vals$stock <- 'global'
      
      #if stocks is null, this will also mean that means/certainty is missing a stocks column
      mean_metric <- as.data.frame(mean_metric)
      mean_metric$stock <- 'global'
      certainty_metric <- as.data.frame(certainty_metric)
      certainty_metric$stock <- 'global'
    }
  } else {
    all_vals <- raw_data
  }
    #PART 1
    # =========================================================================
    # STEP 1: Harmonize Data (Get A, B, C, D counts for both formats)
    # =========================================================================
    if (metric != "sensitivity") {
      # Pivot, bin, and count continuous values
      summary_df <- all_vals %>%
        tidyr::pivot_longer(
          cols = dplyr::any_of(raw_names),
          names_to = "Attribute.Raw",
          values_to = "value"
        ) %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::mutate(
          cat = cut(
            value, 
            breaks = c(-Inf, 1.5, 2.5, 3.5, Inf), 
            labels = c("A", "B", "C", "D")
          )
        ) %>%
        dplyr::group_by(dplyr::across(dplyr::any_of(c("Attribute.Raw", stock_col)))) %>%
        dplyr::summarise(
          A = sum(cat == "A", na.rm = TRUE),
          B = sum(cat == "B", na.rm = TRUE),
          C = sum(cat == "C", na.rm = TRUE),
          D = sum(cat == "D", na.rm = TRUE),
          .groups = "drop"
        )
      
    } else if (metric == "sensitivity") {
      # Sum up existing tally columns
      summary_df <- all_vals %>%
        dplyr::group_by(dplyr::across(dplyr::any_of(c("Attribute.Name", stock_col)))) %>%
        dplyr::summarise(
          A = sum(.data[["Scoring.Rank1"]], na.rm = TRUE),
          B = sum(.data[["Scoring.Rank2"]], na.rm = TRUE),
          C = sum(.data[["Scoring.Rank3"]], na.rm = TRUE),
          D = sum(.data[["Scoring.Rank4"]], na.rm = TRUE),
          .groups = "drop"
        )
    }
    
    # =========================================================================
    # STEP 2: Optional Attribute Mapping (mostly used for continuous)
    # =========================================================================
    if (!is.null(attr_map) && "Attribute.Raw" %in% names(summary_df)) {
      summary_df <- summary_df %>%
        dplyr::mutate(Attribute.Name = attr_map[Attribute.Raw])
    }
    
    # =========================================================================
    # STEP 3: Unified Plot Generation
    # =========================================================================
    # Both formats now have explicit A, B, C, D counts, so we can use geom_col()
    plot_df <- summary_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        p = list(
          data.frame(x = c("A", "B", "C", "D"), y = c(A, B, C, D)) %>%
            ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = x)) +
            ggplot2::geom_col(
              show.legend = FALSE,
              width = 0.95,
              color = 'grey25',
              linewidth = 0.5 
            ) +
            ggplot2::scale_fill_manual(
              values = c("A" = "#008000", "B" = "yellow", "C" = "orange", "D" = "red"),
              drop = FALSE
            ) +
            ggplot2::scale_x_discrete(drop = FALSE) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
            ggplot2::theme_void() +
            ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
        )
      ) %>%
      dplyr::ungroup()
  
  # 4. Pivot Exposure Means and DQ to Long Format
  spMeans_long <- mean_metric %>%
    dplyr::select(
      !!rlang::sym(stock_col),
      dplyr::any_of(raw_names)
    ) %>% # Changed
    tidyr::pivot_longer(
      cols = dplyr::any_of(raw_names), # Changed
      names_to = "Attribute.Raw",
      values_to = "expert.scores"
    ) %>%
    dplyr::mutate(expert.scores = round(as.numeric(expert.scores), 2))
  
  spDQ_long <- certainty_metric %>%
    dplyr::select(
      !!rlang::sym(stock_col),
      dplyr::any_of(raw_names)
    ) %>% # Changed
    tidyr::pivot_longer(
      cols = dplyr::any_of(raw_names), # Changed
      names_to = "Attribute.Raw",
      values_to = "data_quality"
    ) %>%
    dplyr::mutate(data_quality = round(as.numeric(data_quality), 2))
  
  # ---------------------------------------------------------
  # STANDARDIZE STOCK ORDER (Matches short data names to long order list)
  # ---------------------------------------------------------
  raw_stocks <- unique(mean_metric[,stock_col])
  
  # Temporarily translate the short raw names into long names using your stock_key
  long_names_present <- ifelse(
    raw_stocks %in% names(stock_key),
    stock_key[raw_stocks],
    raw_stocks
  )
  
  # Sort the short raw_stocks based on where their long names appear in your stock_order
  ordered_stocks <- raw_stocks[order(
    match(long_names_present, stock_order),
    na.last = TRUE
  )]
  # ---------------------------------------------------------
  
  # 5. Join into Final Summary Table and Pivot WIDER
  tab_wide <- spMeans_long %>%
    dplyr::left_join(spDQ_long, by = c(stock_col, "Attribute.Raw")) %>%
    dplyr::mutate(Attribute.Name = attr_map[Attribute.Raw]) %>%
    dplyr::left_join(
      plot_df,
      by = c("Attribute.Name", stock_col)
    ) %>%
    dplyr::mutate(
      Attribute.Name = factor(Attribute.Name, levels = clean_names)
    ) %>%
    dplyr::select(
      Attribute.Name,
      !!rlang::sym(stock_col),
      expert.scores,
      data_quality,
      p
    ) %>%
    # Lock the column order before pivoting
    dplyr::mutate(
      !!rlang::sym(stock_col) := factor(
        !!rlang::sym(stock_col),
        levels = ordered_stocks
      )
    ) %>% 
    tidyr::pivot_wider(
      names_from = dplyr::all_of(stock_col),
      values_from = c(expert.scores, data_quality, p),
      names_glue = paste0("{.value}_{", stock_col, "}")
    ) %>%
    dplyr::arrange(Attribute.Name)
  
  # 5.5 Extract plot list-columns SAFELY to avoid 'gt' conversion errors
  plot_cols <- grep("^p_", colnames(tab_wide), value = TRUE)
  
  # Only attempt extraction if the columns are STILL lists 
  # (This prevents the code from breaking if you run it twice)
  if (all(sapply(tab_wide[plot_cols], is.list))) {
    
    # Create an empty blank plot for any missing data
    empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()
    
    plot_list_wide <- lapply(tab_wide[plot_cols], function(lst) {
      # Loop through the list and replace NAs/NULLs with an empty plot
      lapply(lst, function(p) {
        if (is.null(p) || (is.logical(p) && is.na(p))) empty_plot else p
      })
    })
  }
  
  # Mutate the dataframe safely to replace plot columns with empty strings
  tab_wide <- tab_wide %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(plot_cols), ~ "")
    )
  
  # 5.6 Build and Append the Summary Row for Sensitivity Data 
  if(metric == 'sensitivity'){
    tab_wide <- tab_wide %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    
    summary_row <- data.frame(Attribute.Name = "Total Sensitivity | Certainty")
    
    for (s in ordered_stocks) {
      # Use %in% instead of == to safely match NA values if a species has no stocks
      s_sens <- mean_metric[[total_sens_col]][mean_metric[[stock_col]] %in% s][1]
      s_cert <- mean_metric[[certainty_col]][mean_metric[[stock_col]] %in% s][1]
      
      sens_word <- switch(
        as.character(s_sens),
        '1' = 'Low',
        '2' = 'Moderate',
        '3' = "High",
        '4' = "Very High",
        'Unknown'
      )
      
      # A zero-width wrapper prevents the 45px column from stretching, 
      # while the inner div spans the full 140px across the 3 columns.
      summary_text <- paste0(
        "<div style='width: 0px; overflow: visible;'>",
        "<div style='width: 140px; text-align: center; white-space: nowrap;'>",
        sens_word, " (", s_sens, ") | ", round(s_cert, 2),
        "</div>",
        "</div>"
      )
      
      # Place it in the FIRST column (expert.scores) so it spans to the right natively
      summary_row[[paste0("expert.scores_", s)]] <- summary_text
      summary_row[[paste0("data_quality_", s)]] <- ""
      summary_row[[paste0("p_", s)]] <- ""
    }
    
    tab_wide <- dplyr::bind_rows(tab_wide, summary_row)
  }
  
  # ---------------------------------------------------------
  # FORCE EXACT COLUMN ORDER (North to South, Grouped by Stock)
  # ---------------------------------------------------------
  desired_cols <- "Attribute.Name"
  for (k in ordered_stocks) {
    desired_cols <- c(
      desired_cols,
      paste0("expert.scores_", k),
      paste0("data_quality_", k),
      paste0("p_", k)
    )
  }
  
  # Reorder the dataframe columns to match the exact custom order
  tab_wide <- tab_wide %>% dplyr::select(dplyr::all_of(desired_cols))
  # ---------------------------------------------------------
  
  # 6. Initialize gt table
  summary.table <- gt::gt(tab_wide) %>%
    gt::tab_header(title = stringr::str_to_title(metric)) %>%
    gt::opt_row_striping() %>%
    gt::cols_align(align = "center", columns = gt::everything()) %>%
    gt::cols_align(align = "left", columns = c("Attribute.Name"))
  
  # Dynamically add spanners for each stock present
  stock_nms <- unique(mean_metric$stock)
  for (k in ordered_stocks) {
    # If a species has no stocks (NA or empty string), skip the spanner entirely
    if (is.na(k) || k == "" || k == "None") {
      next
    }
    
    # Translate the abbreviation if it exists in the map; otherwise use the original string
    display_name <- ifelse(k %in% names(stock_key), stock_key[[k]], k)
    
    wrapped_label <- stringr::str_replace_all(
      stringr::str_wrap(display_name, width = 15),
      pattern = "\n",
      replacement = "  \n"
    )
    
    summary.table <- summary.table %>%
      gt::tab_spanner(
        label = gt::md(wrapped_label),
        columns = gt::ends_with(as.character(k))
      )
  }
  
  if (metric == 'sensitivity') {
    nm <- 'Sensitivity Attributes'
  } else if (metric == 'exposure') {
    nm <- "Exposure Factors"
  } else {
    nm <- "Overall Vulnerability Rank"
  }
  
  # 7. Apply styling, column renaming, and rendering
  summary.table <- summary.table %>%
    # Clean up the Attribute Name header
    gt::cols_label(
      Attribute.Name = nm
    ) %>%
    # Rename grouped columns
    gt::cols_label_with(
      columns = gt::starts_with("expert.scores"),
      fn = ~ ifelse(metric== 'sensitivity', 'Score', "Mean")
    ) %>%
    gt::cols_label_with(
      columns = gt::starts_with("data_quality"),
      fn = ~ ifelse(metric== 'sensitivity', 'Data Quality', "SD")
    ) %>%
    gt::cols_label_with(
      columns = gt::starts_with("p_"),
      fn = ~"Tally"
    ) %>%
    # Force text wrapping by constraining column widths
    gt::cols_width(
      Attribute.Name ~ gt::px(220),
      gt::starts_with("expert.scores") ~ gt::px(55),
      gt::starts_with("data_quality") ~ gt::px(55),
      gt::starts_with("p_") ~ gt::px(55) # Reduced from 90
    ) %>%
    # Make all column headers and spanners bold
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = list(
        gt::cells_column_labels(),
        gt::cells_column_spanners()
      )
    ) %>%
    # Add thin bottom border to stock spanners
    gt::tab_style(
      style = gt::cell_borders(
        sides = "bottom",
        color = "black",
        weight = gt::px(1)
      ),
      locations = gt::cells_column_spanners()
    ) %>%
    # Add thick vertical borders to separate stocks
    gt::tab_style(
      style = gt::cell_borders(
        sides = "left",
        color = "black",
        weight = gt::px(2)
      ),
      locations = gt::cells_body(columns = gt::starts_with("expert.scores"))
    )
  # STOP HERE - NO PIPE (%>%) BEFORE THE FOR LOOP
  
  # Calculate exactly how many rows should receive plots
  num_plot_rows <- length(plot_list_wide[[plot_cols[1]]])
  
  # Inject plots column by column
  for (p_col in plot_cols) {
    summary.table <- local({
      col_name <- p_col
      summary.table %>%
        gt::text_transform(
          locations = gt::cells_body(
            columns = dplyr::all_of(col_name),
            rows = 1:num_plot_rows # Restrict plotting to original rows
          ),
          fn = function(x) {
            purrr::map(plot_list_wide[[col_name]], function(p) {
              # Failsafe: If the plot is NA or missing, return a blank space instead of crashing
              if (is.null(p) || (is.logical(p) && is.na(p[1]))) {
                return("")
              }
              
              gt::ggplot_image(p, height = gt::px(15), aspect_ratio = 3)
            })
          }
        )
    })
  }
  
  # Resume the pipe chain for final table options
  summary.table <- summary.table %>%
    # Parse the custom div strings as HTML instead of raw text
    gt::fmt_markdown(columns = gt::starts_with("expert.scores")) %>%
    
    # 1. Horizontal borders for headers/spanners (Top and Bottom only to prevent doubling)
    gt::tab_style(
      style = gt::cell_borders(
        sides = c("top", "bottom"),
        color = "black",
        weight = gt::px(1)
      ),
      locations = list(
        gt::cells_column_spanners(),
        gt::cells_column_labels()
      )
    ) %>%
    
    # 2. Thick vertical borders (Left side of each stock group) for Body, Headers, & Spanners
    gt::tab_style(
      style = gt::cell_borders(
        sides = "left",
        color = "black",
        weight = gt::px(2)
      ),
      locations = list(
        gt::cells_body(columns = gt::starts_with("expert.scores")),
        gt::cells_column_labels(columns = gt::starts_with("expert.scores")),
        gt::cells_column_spanners() 
      )
    ) %>%
    
    # 3. Thin vertical borders (Right side of inner columns) for Body & Headers
    gt::tab_style(
      style = gt::cell_borders(
        sides = "right",
        color = "black",
        weight = gt::px(1)
      ),
      locations = list(
        gt::cells_body(
          columns = c(gt::starts_with("expert.scores"), gt::starts_with("data_quality")),
          rows = if (metric == "sensitivity") 1:(nrow(tab_wide) - 1) else gt::everything()
        ),
        gt::cells_column_labels(
          columns = c(gt::starts_with("expert.scores"), gt::starts_with("data_quality"))
        )
      )
    ) %>%
    
    # (UPDATED) Thick top line placed above the total section
    gt::tab_style(
      style = gt::cell_borders(
        sides = "top",
        color = "black",
        weight = gt::px(2)
      ),
      locations = gt::cells_body(
        rows = ifelse(metric == 'sensitivity', nrow(tab_wide), nrow(tab_wide) - 1)
      ) 
    )  %>%
    # ... [Keep your other v_align and tab_options here] ...
    
    # 4. Add thick horizontal line to separate the main attributes from the total rows
    gt::tab_style(
      style = gt::cell_borders(
        sides = "top",
        color = "black",
        weight = gt::px(2)
      ),
      locations = gt::cells_body(rows = ifelse(metric == 'sensitivity', nrow(tab_wide), nrow(tab_wide) - 1)) 
    )  %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(rows = nrow(tab_wide)) 
    ) %>%
    # Anchor all cell content to the bottom of the row
    gt::tab_style(
      style = gt::cell_text(v_align = "bottom"),
      locations = gt::cells_body(columns = gt::starts_with("p_"))
    ) %>%
    # Compress padding, center title, and set landscape layout
    gt::tab_options(
      heading.align = "center",
      table.font.size = gt::px(10),
      data_row.padding = gt::px(2),
      heading.padding = gt::px(2),
      column_labels.padding = gt::px(2),
      row.striping.background_color = "#D3D3D3",
      table.border.top.style = "solid",
      table.border.top.width = gt::px(2),
      table.border.top.color = "black",
      table.border.bottom.style = "solid",
      table.border.bottom.width = gt::px(2),
      table.border.bottom.color = "black",
      # APPEND the td overflow rule here:
      table.additional_css = "@page { size: landscape; margin: 0.5in; } td { overflow: visible !important; } img { vertical-align: bottom; margin-bottom: -2px; }"
    ) %>%
    gt::opt_table_lines(extent = "none") %>% gt::opt_table_outline(style = "solid", width = gt::px(2), color = "black")
  
  if(metric == 'exposure'){
    #create flag for which variables are important
    var_imp_nms <- names(var_imp)[var_imp >= imp_threshold]
    
    # 1. Filter the exact list that aligns with your TRUE/FALSE vector
    important_raw_names <- raw_names %in% var_imp_nms
    
    # 2. Translate those specific raw names into their clean display names
    important_clean_names <- attr_map[important_raw_names]
    
    # 3. Append the exact string of your bottom row so it also becomes italic
    important_clean_names <- c(
      important_clean_names,
      clean_names[length(clean_names)] #assumes the string to be italicized is the last clean name 
    )
    
    # Italicize the important attribute names and the important total row
    summary.table <- summary.table %>% 
      gt::tab_style(
        style = gt::cell_text(style = "italic"),
        locations = gt::cells_body(
          rows = Attribute.Name %in% important_clean_names
        )
      ) 
    # (The rogue uppercase/gray styling block was completely removed from here)
  }
  
  if(metric == 'vulnerability'){
    # Italicize the important attribute names and the important total row
    summary.table <- summary.table %>% 
      gt::tab_style(
        style = gt::cell_text(style = "italic"),
        locations = gt::cells_body(
          rows = Attribute.Name %in% clean_names[length(clean_names)] #assumes the string to be italicized is the last clean name
        )
      ) 
  }
  
  if (metric != 'sensitivity') {
    # Both Exposure and Vulnerability now have exactly 2 total rows at the bottom
    bold_rows <- (nrow(tab_wide) - 1):nrow(tab_wide)
    
    summary.table <- summary.table %>%
      # 1. Bold the last two rows entirely
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_body(rows = bold_rows) 
      ) %>%
      # 2. Remove the horizontal line between them to fuse them into one block
      gt::tab_style(
        style = gt::cell_borders(sides = c("top", "bottom"), style = "hidden"),
        locations = gt::cells_body(rows = nrow(tab_wide))
      )
    
  }
  
  # Add a color legend to the bottom of the Vulnerability table
  if (metric == 'vulnerability') {
    
    legend_html <- "
      <div style='text-align: center; padding-top: 5px; font-size: 12px;'>
        <span style='display: inline-block; width: 12px; height: 12px; background-color: green; border: 1px solid black; vertical-align: middle;'></span> <span style='vertical-align: middle;'>Low</span>
        <span style='display: inline-block; width: 12px; height: 12px; background-color: yellow; border: 1px solid black; vertical-align: middle; margin-left: 15px;'></span> <span style='vertical-align: middle;'>Moderate</span>
        <span style='display: inline-block; width: 12px; height: 12px; background-color: orange; border: 1px solid black; vertical-align: middle; margin-left: 15px;'></span> <span style='vertical-align: middle;'>High</span>
        <span style='display: inline-block; width: 12px; height: 12px; background-color: red; border: 1px solid black; vertical-align: middle; margin-left: 15px;'></span> <span style='vertical-align: middle;'>Very High</span>
      </div>
    "
    
    summary.table <- summary.table %>%
      gt::tab_source_note(
        source_note = gt::html(legend_html)
      ) %>%
      # Ensure the source note area integrates seamlessly into the table borders
      gt::tab_options(
        source_notes.background.color = "white",
        source_notes.border.bottom.style = "none",
        source_notes.padding = gt::px(5)
      )
  }
  
  if(save){
    # 8. Save as a cropped, high-res image
    gt::gtsave(
      summary.table,
      paste0(table_dir, '/', metric, '_table_', gsub(' ', '', species), '.png'),
      vwidth = 1500 # Gives it a wide, high-resolution rendering
    )
  }
  
  return(summary.table)
  
}
      
