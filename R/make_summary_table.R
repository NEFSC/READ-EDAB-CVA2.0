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
      ext_vals <- terra::extract(raw_data, stocks)
      ext_vals[[stock_col]] <- stocks$stock_area[ext_vals$ID]
      ext_vals$ID <- NULL
      
      global_vals <- as.data.frame(terra::values(raw_data), check.names = FALSE)
      global_vals[[stock_col]] <- 'global'
      
      all_vals <- dplyr::bind_rows(ext_vals, global_vals)
    } else {
      all_vals <- as.data.frame(terra::values(raw_data), check.names = FALSE)
      all_vals[[stock_col]] <- 'global'
    }
  } else {
    all_vals <- raw_data
  }
  
  # =========================================================================
  # STEP 1.5: Aggressive Formatting to Wide Structure
  # =========================================================================
  
  repair_names <- function(nms, targets) {
    res <- as.character(nms)
    for (i in seq_along(res)) {
      if (res[i] %in% targets) next
      clean_target <- make.names(targets)
      if (res[i] %in% clean_target) {
        res[i] <- targets[clean_target == res[i]][1]
      } else if (res[i] %in% gsub(" ", ".", targets)) {
        res[i] <- targets[gsub(" ", ".", targets) == res[i]][1]
      }
    }
    return(res)
  }
  
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
          df <- tidyr::pivot_wider(df, names_from = dplyr::all_of(name_col), values_from = dplyr::all_of(num_col))
        }
      }
    }
    
    # Inject missing stock columns
    if (!sc %in% colnames(df)) {
      df[[sc]] <- if (!is.null(stocks_obj) && !is.null(stocks_obj$stock_area)) stocks_obj$stock_area[1] else "global"
    }
    return(df)
  }
  
  mean_metric <- force_wide(mean_metric, raw_names, stock_col, stocks)
  certainty_metric <- force_wide(certainty_metric, raw_names, stock_col, stocks)
  colnames(all_vals) <- repair_names(colnames(all_vals), raw_names)
  
  #PART 1
  # =========================================================================
  # STEP 1: Harmonize Data (Get A, B, C, D counts for both formats)
  # =========================================================================
  if (metric != "sensitivity") {
    
    val_cols <- intersect(colnames(all_vals), raw_names)
    if (length(val_cols) == 0) {
      val_cols <- setdiff(colnames(all_vals), c(stock_col, "ID", "x", "y", "geometry", "Species", total_sens_col, certainty_col))
    }
    
    summary_df <- all_vals %>%
      tidyr::pivot_longer(
        cols = dplyr::any_of(val_cols),
        names_to = "Attribute.Raw",
        values_to = "value",
        values_transform = list(value = as.numeric) 
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
  # STEP 2: Optional Attribute Mapping 
  # =========================================================================
  if (!is.null(attr_map) && "Attribute.Raw" %in% names(summary_df)) {
    summary_df <- summary_df %>%
      dplyr::mutate(Attribute.Name = attr_map[Attribute.Raw])
  }
  
  # =========================================================================
  # STEP 3: Unified Plot Generation
  # =========================================================================
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
  
  # =========================================================================
  # STEP 4: Pivot Exposure Means and DQ to Long Format safely
  # =========================================================================
  
  spMeans_long <- mean_metric %>%
    dplyr::select(dplyr::any_of(c(stock_col, raw_names))) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(raw_names),
      names_to = "Attribute.Raw",
      values_to = "expert.scores",
      values_transform = list(expert.scores = as.numeric)
    ) %>%
    dplyr::mutate(expert.scores = round(expert.scores, 2))
  
  spDQ_long <- certainty_metric %>%
    dplyr::select(dplyr::any_of(c(stock_col, raw_names))) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(raw_names), 
      names_to = "Attribute.Raw",
      values_to = "data_quality",
      values_transform = list(data_quality = as.numeric)
    ) %>%
    dplyr::mutate(data_quality = round(data_quality, 2))
  
  # ---------------------------------------------------------
  # STANDARDIZE STOCK ORDER (Matches short data names to long order list)
  # ---------------------------------------------------------
  raw_stocks <- unique(mean_metric[[stock_col]])
  
  long_names_present <- ifelse(
    raw_stocks %in% names(stock_key),
    stock_key[raw_stocks],
    raw_stocks
  )
  
  ordered_stocks <- raw_stocks[order(
    match(long_names_present, stock_order),
    na.last = TRUE
  )]
  # ---------------------------------------------------------
  
  # 5. Join into Final Summary Table and Pivot WIDER
  tab_wide <- spMeans_long %>%
    dplyr::left_join(spDQ_long, by = c(stock_col, "Attribute.Raw")) %>%
    dplyr::mutate(Attribute.Name = unname(ifelse(is.na(attr_map[Attribute.Raw]), Attribute.Raw, attr_map[Attribute.Raw]))) %>%
    dplyr::left_join(
      plot_df,
      by = c("Attribute.Name", stock_col)
    ) %>%
    dplyr::mutate(
      sort_order = match(Attribute.Name, clean_names)
    ) %>%
    dplyr::arrange(sort_order, Attribute.Name) %>%
    dplyr::select(-sort_order) %>%
    dplyr::select(
      Attribute.Name,
      !!rlang::sym(stock_col),
      expert.scores,
      data_quality,
      p
    ) %>%
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
    dplyr::mutate(Attribute.Name = as.character(Attribute.Name))
  
  # =========================================================================
  # STEP 5.1: SUBSET TABLE IF VAR_IMP IS PROVIDED (EXPOSURE ONLY)
  # =========================================================================
  if (!is.null(var_imp) && metric == "exposure") {
    var_imp_nms <- names(var_imp)[var_imp >= imp_threshold]
    
    important_raw_names <- raw_names[raw_names %in% var_imp_nms]
    important_clean_names <- attr_map[important_raw_names]
    
    tab_wide <- tab_wide %>%
      dplyr::filter(Attribute.Name %in% important_clean_names | grepl("Total", Attribute.Name, ignore.case = TRUE))
  }
  
  # =========================================================================
  # STEP 5.2: ENSURE EXPOSURE AND VULNERABILITY ONLY HAVE ONE TOTAL ROW
  # =========================================================================
  if (metric %in% c("exposure", "vulnerability")) {
    total_idx <- grep("Total", tab_wide$Attribute.Name, ignore.case = TRUE)
    if (length(total_idx) > 1) {
      rows_to_drop <- total_idx[-length(total_idx)]
      tab_wide <- tab_wide[-rows_to_drop, ]
    }
  }
  
  # 5.5 Extract plot list-columns SAFELY to avoid 'gt' conversion errors
  plot_cols <- grep("^p_", colnames(tab_wide), value = TRUE)
  
  if (all(sapply(tab_wide[plot_cols], is.list))) {
    empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()
    
    plot_list_wide <- lapply(tab_wide[plot_cols], function(lst) {
      lapply(lst, function(p) {
        if (is.null(p) || (is.logical(p) && is.na(p[1]))) empty_plot else p
      })
    })
  }
  
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
      if (!paste0("expert.scores_", s) %in% colnames(tab_wide)) next
      
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
      
      summary_text <- paste0(
        "<div style='width: 0px; overflow: visible;'>",
        "<div style='width: 175px; text-align: center; white-space: nowrap;'>",
        sens_word, " (", s_sens, ") | ", round(s_cert, 2),
        "</div>",
        "</div>"
      )
      
      summary_row[[paste0("expert.scores_", s)]] <- summary_text
      summary_row[[paste0("data_quality_", s)]] <- ""
      summary_row[[paste0("p_", s)]] <- ""
    }
    
    tab_wide <- dplyr::bind_rows(tab_wide, summary_row)
  }
  
  # 5.7 Add Categorical Text to Total Rows for Exposure and Vulnerability
  if (metric %in% c("exposure", "vulnerability")) {
    total_rows <- grep("Total", tab_wide$Attribute.Name)
    
    for (k in ordered_stocks) {
      exp_col <- paste0("expert.scores_", k)
      
      if (exp_col %in% colnames(tab_wide)) {
        curr_vals <- tab_wide[[exp_col]][total_rows]
        num_vals <- as.numeric(curr_vals)
        
        if (metric == "exposure") {
          cats <- dplyr::case_when(
            is.na(num_vals) ~ "",
            num_vals < 1.5 ~ "Low",
            num_vals < 2.5 ~ "Moderate",
            num_vals < 3.5 ~ "High",
            num_vals >= 3.5 ~ "Very High",
            TRUE ~ ""
          )
        } else if (metric == "vulnerability") {
          cats <- dplyr::case_when(
            is.na(num_vals) ~ "",
            num_vals < 3.5 ~ "Low",
            num_vals < 7.0 ~ "Moderate",
            num_vals < 10.5 ~ "High",
            num_vals >= 10.5 ~ "Very High", 
            TRUE ~ ""
          )
        }
        
        tab_wide[[exp_col]][total_rows] <- ifelse(
          is.na(num_vals),
          curr_vals,
          paste0(cats, "<br>(", curr_vals, ")")
        )
      }
    }
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
  
  tab_wide <- tab_wide %>% dplyr::select(dplyr::any_of(desired_cols))
  
  # 6. Initialize gt table
  summary.table <- gt::gt(tab_wide) %>%
    gt::tab_header(title = stringr::str_to_title(metric)) %>%
    gt::opt_row_striping() %>%
    gt::cols_align(align = "center", columns = gt::everything()) %>%
    gt::cols_align(align = "left", columns = c("Attribute.Name"))
  
  stock_nms <- unique(mean_metric[[stock_col]])
  for (k in ordered_stocks) {
    if (is.na(k) || k == "" || k == "None") {
      next
    }
    
    if (length(ordered_stocks) == 1 && metric %in% c("exposure", "vulnerability")) {
      next
    }
    
    if (!paste0("expert.scores_", k) %in% colnames(tab_wide)) {
      next
    }
    
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
    nm <- "Vulnerability"
  }
  
  # 7. Apply styling, column renaming, and rendering
  summary.table <- summary.table %>%
    gt::cols_label(
      Attribute.Name = nm
    ) %>%
    gt::cols_label_with(
      columns = gt::starts_with("expert.scores"),
      fn = ~ ifelse(metric== 'sensitivity', 'Score', "Mean")
    ) %>%
    gt::cols_label_with(
      columns = gt::starts_with("data_quality"),
      fn = ~ ifelse(metric== 'sensitivity', 'LH DQ', "SD")
    ) %>%
    gt::cols_label_with(
      columns = gt::starts_with("p_"),
      fn = ~"Tally"
    ) %>%
    gt::cols_width(
      Attribute.Name ~ gt::px(220),
      gt::starts_with("expert.scores") ~ gt::px(65),
      gt::starts_with("data_quality") ~ gt::px(55),
      gt::starts_with("p_") ~ gt::px(55) 
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = list(
        gt::cells_column_labels(),
        gt::cells_column_spanners()
      )
    ) %>%
    gt::tab_style(
      style = gt::cell_borders(
        sides = "bottom",
        color = "black",
        weight = gt::px(1)
      ),
      locations = gt::cells_column_spanners()
    ) %>%
    gt::tab_style(
      style = gt::cell_borders(
        sides = "left",
        color = "black",
        weight = gt::px(2)
      ),
      locations = gt::cells_body(columns = gt::starts_with("expert.scores"))
    )
  
  num_plot_rows <- length(plot_list_wide[[plot_cols[1]]])
  
  if (num_plot_rows > 0) {
    for (p_col in plot_cols) {
      summary.table <- local({
        col_name <- p_col
        summary.table %>%
          gt::text_transform(
            locations = gt::cells_body(
              columns = dplyr::all_of(col_name),
              rows = seq_len(num_plot_rows) 
            ),
            fn = function(x) {
              purrr::map(plot_list_wide[[col_name]], function(p) {
                if (is.null(p) || (is.logical(p) && is.na(p[1]))) {
                  return("")
                }
                gt::ggplot_image(p, height = gt::px(15), aspect_ratio = 3)
              })
            }
          )
      })
    }
  }
  
  summary.table <- summary.table %>%
    gt::fmt_markdown(columns = gt::starts_with("expert.scores")) %>%
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
    gt::tab_style(
      style = gt::cell_borders(
        sides = "right",
        color = "black",
        weight = gt::px(1)
      ),
      locations = list(
        gt::cells_body(
          columns = c(gt::starts_with("expert.scores"), gt::starts_with("data_quality")),
          rows = if (metric == "sensitivity") { if (nrow(tab_wide) > 1) seq_len(nrow(tab_wide) - 1) else integer(0) } else gt::everything()
        ),
        gt::cells_column_labels(
          columns = c(gt::starts_with("expert.scores"), gt::starts_with("data_quality"))
        )
      )
    ) %>%
    gt::tab_style(
      style = gt::cell_borders(
        sides = "top",
        color = "black",
        weight = gt::px(2)
      ),
      locations = gt::cells_body(
        rows = if (nrow(tab_wide) > 0) nrow(tab_wide) else integer(0)
      ) 
    )  %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        rows = if (nrow(tab_wide) > 0) nrow(tab_wide) else integer(0)
      ) 
    ) %>%
    gt::tab_style(
      style = gt::cell_text(v_align = "bottom"),
      locations = gt::cells_body(columns = gt::starts_with("p_"))
    ) %>%
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
      table.additional_css = "@page { size: landscape; margin: 0.5in; } td { overflow: visible !important; } img { vertical-align: bottom; margin-bottom: -2px; }"
    ) %>%
    gt::opt_table_lines(extent = "none") %>% gt::opt_table_outline(style = "solid", width = gt::px(2), color = "black")
  
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
      gt::tab_options(
        source_notes.background.color = "white",
        source_notes.border.bottom.style = "none",
        source_notes.padding = gt::px(5)
      )
  }
  
  if(save){
    if(metric == 'sensitivity'){
      gt::gtsave(
        summary.table,
        paste0(table_dir, '/', metric, '_table_', gsub(' ', '', species), '.png'),
        vwidth = 1500 
      )
    } else {
      if(is.null(var_imp)){
        gt::gtsave(
          summary.table,
          paste0(table_dir, '/', metric, '_table_allvar_', gsub(' ', '', species), '.png'),
          vwidth = 1500 
        )
      } else {
        gt::gtsave(
          summary.table,
          paste0(table_dir, '/', metric, '_table_impvar_', gsub(' ', '', species), '.png'),
          vwidth = 1500 
        )
      }
    }
  }
  
  return(summary.table)
}