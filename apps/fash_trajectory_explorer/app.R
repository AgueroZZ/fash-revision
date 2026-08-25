options(stringsAsFactors = FALSE)

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Could not locate the workflowr project root.")
    current <- parent
  }
}

load_exact_object <- function(path, object_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, object_name)) {
    stop("Expected only object ", object_name, " in ", path)
  }
  environment[[object_name]]
}

invisible(lapply(
  c("shiny", "DT", "bslib", "ggplot2", "data.table", "fashr"),
  function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("The ", package, " package is required to run the explorer.")
    }
  }
))

workflowr_root <- find_workflowr_root()
app_directory <- file.path(workflowr_root, "apps", "fash_trajectory_explorer")
analysis_directory <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "fash_linear_real_data_ablation"
)
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(analysis_directory, "fash_linear_real_data_helpers.R"))
source(file.path(app_directory, "explorer_helpers.R"))

explorer_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_trajectory_explorer_mixture_predstep1_penalty10",
  "explorer_index.rds"
)
iwp_fit_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
linear_fit_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_linear_real_data_ablation_mixture_predstep1_penalty10",
  "linear_fit_bf.rds"
)
linear_run_status_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_linear_real_data_ablation_mixture_predstep1_penalty10",
  "run_status.rds"
)
if (!all(file.exists(c(
  explorer_cache_path,
  iwp_fit_path,
  linear_fit_path,
  linear_run_status_path
)))) {
  stop(
    "Explorer inputs are missing. See apps/fash_trajectory_explorer/README.md ",
    "for the cache-build command."
  )
}

linear_run_status <- readRDS(linear_run_status_path)
if (!identical(linear_run_status$status, "complete") ||
    !identical(
      unname(tools::md5sum(linear_fit_path)),
      unname(linear_run_status$fit_md5[["linear_bf"]])
    )) {
  stop("The current-PC mixture run is incomplete or its BF fit changed.")
}

explorer_cache <- readRDS(explorer_cache_path)
index <- explorer_cache$index
presets <- explorer_cache$presets
expected_columns <- c(
  "original_index", "key", "gene_symbol", "gene_id", "variant_id",
  "discovery_status_bf", "iwp1_lfdr_raw", "iwp1_lfdr_bf",
  "linear_lfdr_raw", "linear_lfdr_bf", "strober_linear_pvalue",
  "strober_nonlinear_pvalue"
)
if (!identical(
      explorer_cache$configuration$analysis_id,
      "fash_trajectory_explorer_mixture_predstep1_penalty10"
    ) ||
    !identical(explorer_cache$configuration$pair_count, 1009173L) ||
    !identical(
      explorer_cache$configuration$linear_output_id,
      "fash_linear_real_data_ablation_mixture_predstep1_penalty10"
    ) ||
    !identical(explorer_cache$configuration$linear_prior_mode, "mixture_grid") ||
    !isTRUE(all.equal(
      explorer_cache$configuration$linear_grid,
      default_revision_grid(),
      tolerance = 0
    )) ||
    !isTRUE(all.equal(explorer_cache$configuration$pred_step, 1, tolerance = 0)) ||
    !identical(as.integer(explorer_cache$configuration$penalty), 10L) ||
    nrow(index) != 1009173L || !all(expected_columns %in% names(index)) ||
    nrow(presets) < 1L) {
  stop("The retained explorer index failed structural validation.")
}

# ---------------------------------------------------------------------------
# Startup derivations. Every "significant" decision in this app comes from the
# full-data cumulative-lfdr FDR-0.05 calls stored in the cache, or from the
# reported Strober eFDR at 0.05. Filtering and one-variant-per-gene thinning are
# display operations and never recalibrate FDR inside a subset.
# ---------------------------------------------------------------------------

ALPHA <- explorer_cache$configuration$alpha
if (!isTRUE(all.equal(ALPHA, 0.05))) {
  stop("The explorer cache was built at an unexpected discovery threshold.")
}
METHOD_IDS <- explorer_metric_catalog()$id
METHOD_LABELS <- stats::setNames(explorer_metric_catalog()$label, METHOD_IDS)
METHOD_COLORS <- c(
  iwp1 = "#d7301f",
  linear = "#2c7fb8",
  strober_linear = "#6a51a3",
  strober_nonlinear = "#9e6ebd"
)
TREND_LABEL <- "Weighted linear trend (approx.)"
SERIES_COLORS <- c(
  "IWP1 FASH" = unname(METHOD_COLORS[["iwp1"]]),
  "FASH-linear" = unname(METHOD_COLORS[["linear"]])
)
SERIES_COLORS[TREND_LABEL] <- unname(METHOD_COLORS[["strober_linear"]])
# Hard-wrapped: ggplot2 captions do not wrap, and this one has to fit inside the
# on-screen panel as well as the exported figure.
TREND_CAPTION <- paste(
  "Dashed line: inverse-variance weighted least-squares trend through the",
  "beta-hat estimates.\nStrober et al. test genotype-level data, so this is a",
  "visual approximation, not their fit."
)

index$strober_linear_called <- index$strober_linear_efdr <= ALPHA
index$strober_nonlinear_called <- index$strober_nonlinear_efdr <= ALPHA
index$discovery_status_raw <- make_discovery_status(
  index$iwp1_called_raw,
  index$linear_called_raw
)
if (!identical(
  make_discovery_status(index$iwp1_called_bf, index$linear_called_bf),
  index$discovery_status_bf
)) {
  stop("Recomputed BF discovery status disagrees with the cached column.")
}

call_vectors <- lapply(c(bf = "bf", raw = "raw"), function(version) {
  stats::setNames(
    lapply(METHOD_IDS, function(id) {
      as.logical(index[[explorer_metric_field(id, version, "call")]])
    }),
    METHOD_IDS
  )
})
call_codes <- lapply(call_vectors, make_call_code)
search_fields <- make_explorer_search_fields(index)

full_summary <- summarize_pair_view(index$gene_id, call_vectors$bf)
message(
  "Explorer index loaded: ",
  format(full_summary$pairs, big.mark = ","), " pairs across ",
  format(full_summary$genes, big.mark = ","), " genes."
)

message("Loading retained BF-adjusted IWP1 and compact FASH-linear fits.")
iwp_fit <- load_exact_object(iwp_fit_path, "fash_fit1_update")
linear_fit <- readRDS(linear_fit_path)
validate_compact_linear_mixture_fash(
  linear_fit,
  expected_grid = default_revision_grid(),
  expected_pred_step = 1,
  expected_penalty = 10L
)
if (!identical(names(iwp_fit$fash_data$data_list), index$key) ||
    !isTRUE(linear_fit$bf_adjusted) ||
    !identical(linear_fit$unit_ids, index$key) ||
    !isTRUE(all.equal(
      as.numeric(linear_fit$lfdr),
      index$linear_lfdr_bf,
      tolerance = 0
    ))) {
  stop("The trajectory fits are not aligned with the explorer index.")
}

# One lead-variant-per-gene ordering per metric column, computed on the full
# index the first time it is requested and reused afterwards.
lead_variant_cache <- new.env(parent = emptyenv())
lead_variant_indices <- function(column) {
  if (!exists(column, envir = lead_variant_cache, inherits = FALSE)) {
    selected <- select_lead_variant_per_gene(
      index$gene_id,
      index[[column]],
      index$key
    )
    if (length(selected) != full_summary$genes) {
      stop("Lead-variant selection did not return exactly one pair per gene.")
    }
    assign(column, selected, envir = lead_variant_cache)
  }
  get(column, envir = lead_variant_cache, inherits = FALSE)
}

# Narrow display frames for the browser table, one per lfdr version, built once.
browser_display_cache <- new.env(parent = emptyenv())
browser_display_frame <- function(version) {
  if (!exists(version, envir = browser_display_cache, inherits = FALSE)) {
    columns <- vapply(
      METHOD_IDS, explorer_metric_field, character(1),
      version = version, what = "column"
    )
    frame <- data.frame(
      gene_symbol = index$gene_symbol,
      variant_id = index$variant_id,
      gene_id = index$gene_id,
      status = if (version == "bf") {
        index$discovery_status_bf
      } else {
        index$discovery_status_raw
      },
      iwp1 = index[[columns[["iwp1"]]]],
      linear = index[[columns[["linear"]]]],
      strober_linear = index[[columns[["strober_linear"]]]],
      strober_nonlinear = index[[columns[["strober_nonlinear"]]]],
      calls = call_codes[[version]],
      stringsAsFactors = FALSE
    )
    assign(version, frame, envir = browser_display_cache)
  }
  get(version, envir = browser_display_cache, inherits = FALSE)
}

preset_choices <- lapply(
  split(presets, factor(presets$group, levels = unique(presets$group))),
  function(group) {
    stats::setNames(
      as.character(group$original_index),
      sub("^[^#]*#", "#", group$label)
    )
  }
)
initial_index <- presets$original_index[1L]

significance_choices <- c(
  "Any" = "any",
  "Significant" = "yes",
  "Not significant" = "no"
)
rank_choices <- stats::setNames(METHOD_IDS, paste(
  METHOD_LABELS,
  c("lfdr", "lfdr", "p-value", "p-value")
))
threshold_inputs <- stats::setNames(paste0("max_", METHOD_IDS), METHOD_IDS)

format_count <- function(value) format(as.integer(value), big.mark = ",")

# ---------------------------------------------------------------------------
# User interface.
# ---------------------------------------------------------------------------

significance_control <- function(id) {
  shiny::selectInput(
    paste0("sig_", id),
    METHOD_LABELS[[id]],
    choices = significance_choices,
    selected = "any",
    width = "100%"
  )
}

threshold_control <- function(id) {
  shiny::numericInput(
    threshold_inputs[[id]],
    paste0(
      METHOD_LABELS[[id]], " ",
      if (id %in% c("iwp1", "linear")) "lfdr" else "p", " at most"
    ),
    value = 1,
    min = 0,
    max = 1,
    step = 0.01,
    width = "100%"
  )
}

trajectory_tab <- bslib::nav_panel(
  title = "Trajectory explorer",
  value = "explorer",
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 350,
      shiny::tags$p(
        class = "explorer-intro",
        paste0(
          "Search all ", format_count(full_summary$pairs),
          " tested pairs by HGNC symbol, Ensembl ID, or rsID. Two terms are",
          " combined, for example \"GPR78 rs4583742\"."
        )
      ),
      shiny::textInput(
        "query",
        "Gene or variant",
        placeholder = "e.g. GPR78, ENSG00000155269, rs4583742"
      ),
      shiny::div(
        class = "search-row",
        shiny::actionButton("search", "Search", class = "btn-primary"),
        shiny::actionButton("previous_pair", "◀", class = "btn-step", title = "Previous result"),
        shiny::actionButton("next_pair", "▶", class = "btn-step", title = "Next result")
      ),
      shiny::uiOutput("search_status"),
      shiny::tags$hr(),
      shiny::selectInput(
        "preset",
        "Curated examples",
        choices = preset_choices,
        selected = as.character(initial_index),
        width = "100%"
      ),
      shiny::checkboxGroupInput(
        "layers",
        "Plot layers",
        choices = c(
          "Observed estimates (plus or minus 2 SE)" = "observed",
          "IWP1 FASH posterior" = "iwp1",
          "FASH-linear posterior" = "linear",
          "95% posterior intervals" = "intervals",
          "Weighted linear trend (Strober-style approx.)" = "trend"
        ),
        selected = c("observed", "iwp1", "linear", "intervals")
      ),
      shiny::tags$div(
        class = "explorer-note explorer-note-tight",
        shiny::tags$strong("The trend layer is an approximation,"),
        " not Strober's genotype-level fit. The exported figure repeats this.",
        " First selection of a pair takes a few seconds, then it is cached."
      ),
      shiny::tags$hr(),
      shiny::downloadButton("download_plot", "Figure (PNG)", class = "btn-light w-100"),
      shiny::downloadButton("download_pair", "This pair (CSV)", class = "btn-light w-100 mt-2")
    ),
    # A plain wrapper keeps layout_sidebar from stretching the result card to
    # fill the panel height, which left a large empty block under the table.
    shiny::div(
      class = "tab-body",
      bslib::card(
        bslib::card_header("Search results — select one row to update the plot"),
        # fill = FALSE on both the card and the output: inside a fill container
        # DT switches on fillContainer and pads the table body with dead space.
        DT::DTOutput("search_results", fill = FALSE),
        fill = FALSE,
        class = "results-card"
      ),
      bslib::layout_columns(
        col_widths = c(8, 4),
        bslib::card(
          full_screen = TRUE,
          bslib::card_header(shiny::uiOutput("selected_pair_header")),
          shiny::plotOutput("trajectory", height = "520px")
        ),
        bslib::card(
          bslib::card_header("Matched evidence"),
          shiny::uiOutput("evidence_cards"),
          shiny::tags$div(
            class = "metric-note",
            "Large number is the headline statistic; small number is its",
            " companion (raw lfdr for FASH, eFDR for Strober). FASH calls use the",
            " cumulative-lfdr FDR 0.05 rule on all",
            paste0(" ", format_count(full_summary$pairs), " pairs;"),
            " Strober calls use the reported eFDR at 0.05. Posterior curves are",
            " always the BF-adjusted fits."
          )
        )
      )
    )
  )
)

browser_tab <- bslib::nav_panel(
  title = "All-pair browser",
  value = "browser",
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 350,
      shiny::radioButtons(
        "browser_version",
        "FASH lfdr version",
        choices = c("BF-adjusted" = "bf", "Raw" = "raw"),
        selected = "bf",
        inline = TRUE
      ),
      shiny::selectInput(
        "rank_by",
        "Rank by",
        choices = rank_choices,
        selected = "iwp1",
        width = "100%"
      ),
      shiny::textInput(
        "browser_query",
        "Gene or variant contains",
        placeholder = "symbol, Ensembl ID, or rsID"
      ),
      shiny::tags$hr(),
      shiny::checkboxInput(
        "one_per_gene",
        "Keep one variant per gene",
        value = FALSE
      ),
      shiny::conditionalPanel(
        condition = "input.one_per_gene == true",
        shiny::selectInput(
          "lead_metric",
          "Lead variant chosen by",
          choices = c("Same as \"Rank by\"" = "rank", rank_choices),
          selected = "rank",
          width = "100%"
        ),
        shiny::tags$div(
          class = "explorer-note explorer-note-tight",
          "Thinning runs on the full index before any filter, so a gene's lead",
          " variant is its most significant variant among all tested variants.",
          " Ties break on ascending pair key."
        )
      ),
      shiny::tags$hr(),
      shiny::tags$div(class = "sidebar-label", "Significance filters"),
      lapply(METHOD_IDS, significance_control),
      shiny::selectInput(
        "discordance",
        "Disagreement",
        choices = c(
          "Any" = "any",
          "IWP1 vs FASH-linear disagree" = "fash_models",
          "Strober linear vs nonlinear disagree" = "strober_tests",
          "FASH (either) vs Strober (either) disagree" = "fash_vs_strober"
        ),
        selected = "any",
        width = "100%"
      ),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Statistic thresholds",
          lapply(METHOD_IDS, threshold_control)
        )
      ),
      shiny::actionButton("reset_filters", "Reset filters", class = "btn-light w-100"),
      shiny::tags$div(
        class = "explorer-note",
        shiny::tags$strong("Filters never recalibrate FDR."),
        " Significance always reflects the full-data calls over all",
        paste0(" ", format_count(full_summary$pairs), " pairs."),
        " Counts below describe the current view only."
      )
    ),
    shiny::div(
      class = "tab-body",
      shiny::uiOutput("view_summary"),
      bslib::layout_columns(
      col_widths = c(8, 4),
      bslib::card(
        bslib::card_header(
          shiny::uiOutput("browser_table_header"),
          class = "d-flex justify-content-between align-items-center"
        ),
        DT::DTOutput("browser_table", fill = FALSE),
        shiny::tags$div(
          class = "metric-note metric-note-tight",
          shiny::span(class = "call-dots-legend", "Calls column, left to right:"),
          shiny::span(class = "call-dot call-dot-iwp1 is-on"), " IWP1 FASH  ",
          shiny::span(class = "call-dot call-dot-linear is-on"), " FASH-linear  ",
          shiny::span(class = "call-dot call-dot-strober-linear is-on"), " Strober linear  ",
          shiny::span(class = "call-dot call-dot-strober-nonlinear is-on"), " Strober nonlinear.",
          " Filled means significant."
        ),
        shiny::downloadButton("download_view", "Download current view (CSV)", class = "btn-light")
      ),
      bslib::card(
        bslib::card_header("Selected pair"),
        shiny::uiOutput("browser_selection_header"),
        shiny::plotOutput("browser_trajectory", height = "300px"),
        shiny::actionButton(
          "open_in_explorer",
          "Open in trajectory explorer",
          class = "btn-primary w-100"
        ),
        bslib::card_body(
          class = "concordance-body",
          shiny::tags$div(
            class = "sidebar-label",
            "Co-significant pairs in view"
          ),
          shiny::tags$div(
            class = "concordance-hint",
            "Off-diagonal cells count pairs called by both methods; the",
            " diagonal is each method's total in the current view."
          ),
          shiny::uiOutput("concordance")
        )
      )
      )
    )
  )
)

ui <- bslib::page_navbar(
  id = "main_nav",
  title = shiny::div(
    class = "explorer-title",
    "FASH trajectory explorer",
    shiny::span("Time-specific-PC real-data comparison")
  ),
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  header = shiny::tags$head(
    shiny::tags$link(rel = "stylesheet", type = "text/css", href = "explorer.css")
  ),
  trajectory_tab,
  browser_tab
)

# ---------------------------------------------------------------------------
# Plot construction shared by both tabs.
# ---------------------------------------------------------------------------

build_trajectory_plot <- function(row,
                                  trajectory,
                                  layers,
                                  base_size = 13,
                                  compact = FALSE,
                                  annotate = TRUE) {
  plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.45)

  if ("intervals" %in% layers && "iwp1" %in% layers) {
    plot <- plot + ggplot2::geom_ribbon(
      data = trajectory$iwp1,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = method),
      alpha = 0.14
    )
  }
  if ("intervals" %in% layers && "linear" %in% layers) {
    plot <- plot + ggplot2::geom_ribbon(
      data = trajectory$linear,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = method),
      alpha = 0.14
    )
  }
  if ("iwp1" %in% layers) {
    plot <- plot + ggplot2::geom_line(
      data = trajectory$iwp1,
      ggplot2::aes(x = time, y = posterior_mean, color = method, linetype = method),
      linewidth = 1.15
    )
  }
  if ("linear" %in% layers) {
    plot <- plot + ggplot2::geom_line(
      data = trajectory$linear,
      ggplot2::aes(x = time, y = posterior_mean, color = method, linetype = method),
      linewidth = 1.15
    )
  }
  if ("trend" %in% layers) {
    plot <- plot + ggplot2::geom_line(
      data = trajectory$trend,
      ggplot2::aes(x = time, y = fitted, color = method, linetype = method),
      linewidth = 0.95
    )
  }
  if ("observed" %in% layers) {
    plot <- plot +
      ggplot2::geom_errorbar(
        data = trajectory$observed,
        ggplot2::aes(
          x = time,
          ymin = beta - 2 * standard_error,
          ymax = beta + 2 * standard_error
        ),
        width = 0.12,
        color = "grey30",
        linewidth = 0.45
      ) +
      ggplot2::geom_point(
        data = trajectory$observed,
        ggplot2::aes(x = time, y = beta),
        size = if (compact) 1.4 else 2,
        color = "black"
      )
  }

  line_types <- stats::setNames(
    ifelse(names(SERIES_COLORS) == TREND_LABEL, "22", "solid"),
    names(SERIES_COLORS)
  )
  plot +
    ggplot2::scale_color_manual(values = SERIES_COLORS, name = NULL) +
    ggplot2::scale_fill_manual(values = SERIES_COLORS, name = NULL, guide = "none") +
    ggplot2::scale_linetype_manual(values = line_types, name = NULL) +
    ggplot2::labs(
      title = if (annotate) paste(row$gene_symbol, "/", row$variant_id) else NULL,
      subtitle = if (annotate) make_explorer_plot_subtitle(row) else NULL,
      caption = if ("trend" %in% layers && !compact) TREND_CAPTION else NULL,
      x = "Differentiation time",
      y = "Estimated eQTL effect"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = 0),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(lineheight = 1.15),
      plot.caption = ggplot2::element_text(
        hjust = 0, size = base_size - 4, color = "grey35", lineheight = 1.15
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# The dots are rendered client-side from a 4-bit call code so that a filtered
# view of a million rows never has to materialise a million HTML strings.
call_dots_renderer <- DT::JS(
  "function(data, type, row, meta) {",
  "  if (type !== 'display') { return data; }",
  "  var labels = ['IWP1 FASH', 'FASH-linear', 'Strober linear', 'Strober nonlinear'];",
  "  var keys = ['iwp1', 'linear', 'strober-linear', 'strober-nonlinear'];",
  "  var out = '<span class=\"call-dots\">';",
  "  for (var i = 0; i < 4; i++) {",
  "    var on = (data & (1 << i)) > 0;",
  "    out += '<span class=\"call-dot call-dot-' + keys[i] + (on ? ' is-on' : '') +",
  "      '\" title=\"' + labels[i] + ': ' + (on ? 'significant' : 'not significant') + '\"></span>';",
  "  }",
  "  return out + '</span>';",
  "}"
)

# ---------------------------------------------------------------------------
# Server.
# ---------------------------------------------------------------------------

server <- function(input, output, session) {
  selected_index <- shiny::reactiveVal(initial_index)
  initial_results <- index[unique(presets$original_index), , drop = FALSE]
  search_results <- shiny::reactiveVal(initial_results)
  search_message <- shiny::reactiveVal(
    paste("Showing", nrow(initial_results), "curated examples.")
  )
  trajectory_cache <- new.env(parent = emptyenv())

  # -- Trajectory tab: search -------------------------------------------------

  run_search <- function(query, notify) {
    tokens <- split_explorer_tokens(query)
    if (length(tokens) == 0L) {
      search_results(initial_results)
      search_message(paste("Showing", nrow(initial_results), "curated examples."))
      return(invisible(NULL))
    }
    result <- search_explorer_index(index, query, limit = 200L, fields = search_fields)
    search_results(result)
    if (nrow(result) == 0L) {
      search_message("No matching gene-variant pair.")
      if (notify) {
        shiny::showNotification(
          "No matching gene-variant pair was found.",
          type = "warning"
        )
      }
    } else {
      search_message(paste0(
        "Showing ", nrow(result), if (nrow(result) == 200L) " of the top" else "",
        " ranked ", if (nrow(result) == 1L) "match." else "matches."
      ))
    }
    invisible(NULL)
  }

  debounced_query <- shiny::debounce(shiny::reactive(input$query), 350)
  shiny::observeEvent(debounced_query(), {
    run_search(debounced_query(), notify = FALSE)
  }, ignoreInit = TRUE)
  shiny::observeEvent(input$search, {
    run_search(input$query, notify = TRUE)
  }, ignoreInit = TRUE)

  output$search_status <- shiny::renderUI({
    shiny::tags$div(class = "search-status", search_message())
  })

  shiny::observeEvent(input$preset, {
    selected_index(as.integer(input$preset))
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$search_results_rows_selected, {
    selected_row_number <- input$search_results_rows_selected
    result <- search_results()
    if (length(selected_row_number) == 1L && selected_row_number <= nrow(result)) {
      selected_index(result$original_index[selected_row_number])
    }
  }, ignoreInit = TRUE)

  step_selection <- function(offset) {
    result <- search_results()
    if (nrow(result) == 0L) return(invisible(NULL))
    position <- match(selected_index(), result$original_index)
    position <- if (is.na(position)) 1L else position + offset
    position <- max(1L, min(nrow(result), position))
    selected_index(result$original_index[position])
    DT::selectRows(DT::dataTableProxy("search_results"), position)
    invisible(NULL)
  }
  shiny::observeEvent(input$previous_pair, step_selection(-1L), ignoreInit = TRUE)
  shiny::observeEvent(input$next_pair, step_selection(1L), ignoreInit = TRUE)

  # -- Trajectory tab: selected pair -----------------------------------------

  selected_row <- shiny::reactive({
    index[selected_index(), , drop = FALSE]
  })

  selected_trajectory <- shiny::reactive({
    original_index <- selected_index()
    cache_key <- as.character(original_index)
    if (exists(cache_key, envir = trajectory_cache, inherits = FALSE)) {
      return(get(cache_key, envir = trajectory_cache, inherits = FALSE))
    }

    shiny::withProgress(message = "Computing posterior trajectories", value = 0.2, {
      dataset <- iwp_fit$fash_data$data_list[[original_index]]
      standard_error <- as.numeric(iwp_fit$fash_data$S[[original_index]])
      grid <- seq(0, 15, by = 0.1)
      set.seed(20260810L + original_index %% 1000000L)
      iwp_prediction <- stats::predict(
        iwp_fit,
        index = original_index,
        smooth_var = grid
      )
      shiny::incProgress(0.45)
      linear_prediction <- extract_linear_mixture_posterior_plot_data(
        dataset = dataset,
        standard_error = standard_error,
        fit = linear_fit,
        unit_index = original_index,
        grid = grid,
        sample_size = 10000L,
        seed = 20260810L + original_index %% 1000000L
      )
      trend <- fit_weighted_linear_trend(
        time = as.numeric(dataset$x),
        beta = as.numeric(dataset$y),
        standard_error = standard_error,
        grid = grid
      )
      result <- list(
        observed = data.frame(
          time = as.numeric(dataset$x),
          beta = as.numeric(dataset$y),
          standard_error = standard_error,
          stringsAsFactors = FALSE
        ),
        iwp1 = data.frame(
          method = "IWP1 FASH",
          time = as.numeric(iwp_prediction$x),
          posterior_mean = as.numeric(iwp_prediction$mean),
          lower = as.numeric(iwp_prediction$lower),
          upper = as.numeric(iwp_prediction$upper),
          stringsAsFactors = FALSE
        ),
        linear = transform(linear_prediction, method = "FASH-linear"),
        trend = transform(trend, method = TREND_LABEL)
      )
      assign(cache_key, result, envir = trajectory_cache)
      shiny::incProgress(0.35)
      result
    })
  })

  output$selected_pair_header <- shiny::renderUI({
    row <- selected_row()
    shiny::tags$div(
      shiny::tags$strong(paste(row$gene_symbol, "/", row$variant_id)),
      shiny::tags$span(
        class = paste0(
          "status-badge status-",
          gsub("[^a-z]+", "-", tolower(row$discovery_status_bf))
        ),
        row$discovery_status_bf
      ),
      shiny::tags$small(row$gene_id)
    )
  })

  output$evidence_cards <- shiny::renderUI({
    cards <- make_explorer_metric_cards(selected_row(), version = "bf", alpha = ALPHA)
    shiny::tags$div(
      class = "evidence-grid",
      lapply(seq_len(nrow(cards)), function(position) {
        shiny::tags$div(
          class = paste0(
            "evidence-card evidence-", gsub("_", "-", cards$id[position]),
            if (cards$called[position]) " is-called" else ""
          ),
          shiny::tags$div(class = "evidence-method", cards$label[position]),
          shiny::tags$div(
            class = "evidence-primary",
            format_explorer_number(cards$primary_value[position], 3L)
          ),
          shiny::tags$div(
            class = "evidence-primary-label",
            cards$primary_label[position]
          ),
          shiny::tags$div(
            class = "evidence-secondary",
            paste0(
              cards$secondary_label[position], " ",
              format_explorer_number(cards$secondary_value[position], 3L)
            )
          ),
          shiny::tags$span(
            class = paste0(
              "evidence-chip ",
              if (cards$called[position]) "chip-yes" else "chip-no"
            ),
            if (cards$called[position]) "significant" else "not significant"
          )
        )
      })
    )
  })

  output$search_results <- DT::renderDT({
    result <- search_results()
    display <- result[, c(
      "original_index", "gene_symbol", "gene_id", "variant_id",
      "discovery_status_bf", "iwp1_lfdr_bf", "linear_lfdr_bf",
      "strober_linear_pvalue", "strober_nonlinear_pvalue"
    ), drop = FALSE]
    names(display) <- c(
      "Index", "Symbol", "Ensembl", "Variant", "BF call",
      "IWP1 lfdr", "Linear lfdr", "Strober linear p", "Strober nonlinear p"
    )
    table <- DT::datatable(
      display,
      rownames = FALSE,
      selection = "single",
      filter = "top",
      options = list(
        pageLength = 8,
        lengthMenu = c(8, 25, 50, 100),
        scrollX = TRUE,
        order = list(list(5, "asc"), list(6, "asc")),
        columnDefs = list(list(visible = FALSE, targets = 0))
      )
    )
    DT::formatSignif(table, columns = 6:9, digits = 3)
  })

  # On screen the card header and the evidence cards already name the pair and
  # its statistics, so the figure keeps that space for the data. The downloaded
  # PNG has to stand alone, so it carries the title and subtitle.
  build_current_plot <- function(annotate) {
    build_trajectory_plot(
      row = selected_row(),
      trajectory = selected_trajectory(),
      layers = input$layers,
      annotate = annotate
    )
  }

  output$trajectory <- shiny::renderPlot(
    build_current_plot(annotate = FALSE),
    res = 120
  )

  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      row <- selected_row()
      paste0("fash_trajectory_", row$gene_symbol, "_", row$variant_id, ".png")
    },
    content = function(file) {
      ggplot2::ggsave(
        file,
        plot = build_current_plot(annotate = TRUE),
        width = 8, height = 6, dpi = 200
      )
    }
  )

  output$download_pair <- shiny::downloadHandler(
    filename = function() {
      row <- selected_row()
      paste0("fash_pair_", row$gene_symbol, "_", row$variant_id, ".csv")
    },
    content = function(file) {
      row <- selected_row()
      trajectory <- selected_trajectory()
      curves <- rbind(
        data.frame(
          series = "observed",
          time = trajectory$observed$time,
          value = trajectory$observed$beta,
          lower = trajectory$observed$beta - 2 * trajectory$observed$standard_error,
          upper = trajectory$observed$beta + 2 * trajectory$observed$standard_error,
          stringsAsFactors = FALSE
        ),
        data.frame(
          series = "iwp1_posterior_mean",
          time = trajectory$iwp1$time,
          value = trajectory$iwp1$posterior_mean,
          lower = trajectory$iwp1$lower,
          upper = trajectory$iwp1$upper,
          stringsAsFactors = FALSE
        ),
        data.frame(
          series = "linear_posterior_mean",
          time = trajectory$linear$time,
          value = trajectory$linear$posterior_mean,
          lower = trajectory$linear$lower,
          upper = trajectory$linear$upper,
          stringsAsFactors = FALSE
        ),
        data.frame(
          series = "weighted_linear_trend_approximation",
          time = trajectory$trend$time,
          value = trajectory$trend$fitted,
          lower = NA_real_,
          upper = NA_real_,
          stringsAsFactors = FALSE
        )
      )
      curves$gene_symbol <- row$gene_symbol
      curves$gene_id <- row$gene_id
      curves$variant_id <- row$variant_id
      evidence <- make_explorer_metrics(row, alpha = ALPHA)
      connection <- file(file, open = "wt")
      on.exit(close(connection), add = TRUE)
      writeLines(
        c(
          paste0("# pair: ", row$gene_symbol, " / ", row$variant_id),
          paste0("# key: ", row$key),
          "# evidence (full-data calls at FDR 0.05):",
          paste0(
            "#   ", evidence$Analysis, " ", evidence$Metric, " = ",
            evidence$Value, " | adjusted = ", evidence$`Adjusted metric`,
            " | called = ", evidence$`Called at 0.05`
          ),
          paste0("# ", TREND_CAPTION)
        ),
        connection
      )
      utils::write.csv(curves, connection, row.names = FALSE)
    }
  )

  # -- Browser tab: filtering ------------------------------------------------

  browser_version <- shiny::reactive({
    if (identical(input$browser_version, "raw")) "raw" else "bf"
  })
  rank_metric <- shiny::reactive({
    if (is.null(input$rank_by)) "iwp1" else input$rank_by
  })
  rank_column <- shiny::reactive({
    explorer_metric_field(rank_metric(), browser_version(), "column")
  })
  lead_metric <- shiny::reactive({
    choice <- input$lead_metric
    if (is.null(choice) || identical(choice, "rank")) rank_metric() else choice
  })
  debounced_browser_query <- shiny::debounce(shiny::reactive(input$browser_query), 400)

  shiny::observeEvent(input$reset_filters, {
    shiny::updateTextInput(session, "browser_query", value = "")
    shiny::updateCheckboxInput(session, "one_per_gene", value = FALSE)
    shiny::updateSelectInput(session, "discordance", selected = "any")
    for (id in METHOD_IDS) {
      shiny::updateSelectInput(session, paste0("sig_", id), selected = "any")
      shiny::updateNumericInput(session, threshold_inputs[[id]], value = 1)
    }
  })

  view_indices <- shiny::reactive({
    version <- browser_version()
    calls <- call_vectors[[version]]
    keep <- rep(TRUE, nrow(index))

    if (isTRUE(input$one_per_gene)) {
      lead <- lead_variant_indices(
        explorer_metric_field(lead_metric(), version, "column")
      )
      keep <- logical(nrow(index))
      keep[lead] <- TRUE
    }
    keep <- keep & match_explorer_text(search_fields, debounced_browser_query())
    for (id in METHOD_IDS) {
      mode <- input[[paste0("sig_", id)]]
      if (!is.null(mode) && !identical(mode, "any")) {
        keep <- keep & apply_significance_filter(calls[[id]], mode)
      }
    }
    if (!is.null(input$discordance) && !identical(input$discordance, "any")) {
      keep <- keep & apply_discordance_filter(calls, input$discordance)
    }
    for (id in METHOD_IDS) {
      threshold <- sanitize_threshold(input[[threshold_inputs[[id]]]])
      if (threshold < 1) {
        column <- explorer_metric_field(id, version, "column")
        keep <- keep & (index[[column]] <= threshold)
      }
    }
    which(keep)
  })

  # Radix ordering is stable, so equal statistics keep ascending index order.
  ranked_indices <- shiny::reactive({
    indices <- view_indices()
    if (length(indices) == 0L) return(indices)
    indices[order(index[[rank_column()]][indices], method = "radix")]
  })

  view_summary <- shiny::reactive({
    indices <- ranked_indices()
    summarize_pair_view(
      index$gene_id[indices],
      lapply(call_vectors[[browser_version()]], `[`, indices)
    )
  })

  output$view_summary <- shiny::renderUI({
    summary <- view_summary()
    tiles <- list(
      list(label = "Pairs in view", value = summary$pairs, total = full_summary$pairs, id = "pairs"),
      list(label = "Genes in view", value = summary$genes, total = full_summary$genes, id = "genes")
    )
    for (id in METHOD_IDS) {
      tiles[[length(tiles) + 1L]] <- list(
        label = paste(METHOD_LABELS[[id]], "significant"),
        value = summary$significant[[id]],
        total = NULL,
        id = id
      )
    }
    shiny::tags$div(
      class = "summary-strip",
      lapply(tiles, function(tile) {
        shiny::tags$div(
          class = paste0("summary-tile summary-", gsub("_", "-", tile$id)),
          shiny::tags$div(class = "summary-value", format_count(tile$value)),
          shiny::tags$div(class = "summary-label", tile$label),
          if (!is.null(tile$total)) {
            shiny::tags$div(
              class = "summary-total",
              paste0("of ", format_count(tile$total))
            )
          }
        )
      })
    )
  })

  output$concordance <- shiny::renderUI({
    joint <- view_summary()$joint
    labels <- unname(METHOD_LABELS[rownames(joint)])
    shiny::tags$table(
      class = "concordance-table",
      shiny::tags$thead(
        shiny::tags$tr(
          shiny::tags$th(""),
          lapply(labels, function(label) shiny::tags$th(label))
        )
      ),
      shiny::tags$tbody(
        lapply(seq_along(labels), function(i) {
          shiny::tags$tr(
            shiny::tags$th(labels[i]),
            lapply(seq_along(labels), function(j) {
              if (j > i) {
                shiny::tags$td(class = "concordance-empty", "")
              } else {
                shiny::tags$td(
                  class = if (i == j) "concordance-diagonal" else "concordance-cell",
                  format_count(joint[i, j])
                )
              }
            })
          )
        })
      )
    )
  })

  output$browser_table_header <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$span(
        class = "table-title",
        paste0(
          format_count(length(ranked_indices())), " pairs, ranked by ",
          explorer_metric_label(rank_metric(), browser_version())
        )
      ),
      shiny::tags$span(
        class = "table-hint",
        "Click a row to plot it; click a header to re-sort."
      )
    )
  })

  output$browser_table <- DT::renderDT({
    indices <- ranked_indices()
    version <- browser_version()
    display <- browser_display_frame(version)[indices, , drop = FALSE]
    names(display) <- c(
      "Symbol", "Variant", "Ensembl", "Status",
      explorer_metric_label("iwp1", version),
      explorer_metric_label("linear", version),
      "Strober linear p", "Strober nonlinear p", "Calls"
    )
    table <- DT::datatable(
      display,
      rownames = FALSE,
      selection = "single",
      filter = "none",
      escape = FALSE,
      options = list(
        # No DataTables search box: it would filter the page without updating
        # the summary tiles. The sidebar text filter does both.
        dom = "lrtip",
        pageLength = 15,
        lengthMenu = c(15, 25, 50, 100),
        scrollX = TRUE,
        order = list(),
        deferRender = TRUE,
        columnDefs = list(
          list(targets = 8, render = call_dots_renderer, orderable = FALSE,
               className = "dt-center")
        )
      )
    )
    DT::formatSignif(table, columns = 5:8, digits = 3)
  }, server = TRUE)

  shiny::observeEvent(input$browser_table_rows_selected, {
    position <- input$browser_table_rows_selected
    indices <- ranked_indices()
    if (length(position) == 1L && position <= length(indices)) {
      selected_index(indices[position])
    }
  }, ignoreInit = TRUE)

  output$browser_selection_header <- shiny::renderUI({
    row <- selected_row()
    shiny::tags$div(
      class = "browser-selection",
      shiny::tags$strong(paste(row$gene_symbol, "/", row$variant_id)),
      shiny::tags$span(
        class = paste0(
          "status-badge status-",
          gsub("[^a-z]+", "-", tolower(row$discovery_status_bf))
        ),
        row$discovery_status_bf
      ),
      shiny::tags$div(class = "browser-selection-key", row$key)
    )
  })

  output$browser_trajectory <- shiny::renderPlot({
    build_trajectory_plot(
      row = selected_row(),
      trajectory = selected_trajectory(),
      layers = c("observed", "iwp1", "linear"),
      base_size = 11,
      compact = TRUE,
      annotate = FALSE
    )
  }, res = 100)

  shiny::observeEvent(input$open_in_explorer, {
    bslib::nav_select("main_nav", selected = "explorer", session = session)
  })

  output$download_view <- shiny::downloadHandler(
    filename = function() {
      paste0("fash_pair_view_", length(ranked_indices()), "_pairs.csv")
    },
    content = function(file) {
      indices <- ranked_indices()
      shiny::withProgress(message = "Writing the current view", value = 0.4, {
        export <- index[indices, , drop = FALSE]
        export$call_code <- call_codes[[browser_version()]][indices]
        export$rank_metric <- rank_column()
        data.table::fwrite(export, file)
      })
    }
  )
}

shiny::shinyApp(ui, server)
