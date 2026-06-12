library(shiny)

resolve_app_dir <- function() {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- c(cwd, file.path(cwd, "shiny_app"), file.path(cwd, "apps", "shiny_app"))
  app_dirs <- candidates[file.exists(file.path(candidates, "app.R"))]
  if (!length(app_dirs)) {
    stop("Could not find apps/shiny_app/app.R from the current working directory.")
  }
  normalizePath(app_dirs[[1]], winslash = "/", mustWork = TRUE)
}

count_assets <- function(path, pattern = NULL) {
  if (!dir.exists(path)) {
    return(0L)
  }
  length(list.files(path, pattern = pattern, full.names = FALSE))
}

list_assets <- function(path, label) {
  if (!dir.exists(path)) {
    return(data.frame(section = label, file = character(0), stringsAsFactors = FALSE))
  }
  files <- list.files(path, full.names = FALSE)
  if (!length(files)) {
    return(data.frame(section = label, file = "(empty)", stringsAsFactors = FALSE))
  }
  data.frame(section = label, file = files, stringsAsFactors = FALSE)
}

app_dir <- resolve_app_dir()
repo_root <- normalizePath(file.path(app_dir, "..", ".."), winslash = "/", mustWork = TRUE)

ui <- fluidPage(
  tags$head(
    tags$title("Pertussis Visualization Workspace"),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  div(
    class = "page-shell",
    div(
      class = "hero",
      div(
        class = "hero-copy",
        tags$div(class = "eyebrow", "Pertussis Genomics Workspace"),
        tags$h1("Static figures, manuscript assets, and the interactive roadmap in one place."),
        tags$p(
          "This scaffold reads the repository structure directly so collaborators can see what is ready, what is archived, and what still belongs to future work."
        )
      ),
      div(
        class = "hero-note",
        tags$h2("Current delivery model"),
        tags$p("Static publication figures remain the primary output. The Shiny layer is a navigation and inspection surface until curated datasets are finalized.")
      )
    ),
    uiOutput("summary_cards"),
    div(
      class = "panel-grid",
      div(
        class = "panel-card wide",
        tags$h3("Figure asset inventory"),
        tags$p("The app lists the figure workspace exactly as it exists on disk."),
        tableOutput("asset_table")
      ),
      div(
        class = "panel-card",
        tags$h3("Manuscript extracts"),
        tags$p("Frozen TSV extracts available for manuscript-facing figures."),
        tableOutput("figure_data_table")
      ),
      div(
        class = "panel-card",
        tags$h3("Module status"),
        tags$div(class = "module-state active", tags$strong("Static figures"), tags$span("Active")),
        tags$div(class = "module-state active", tags$strong("Manuscript assets"), tags$span("Active")),
        tags$div(class = "module-state pending", tags$strong("Public-health ingestion"), tags$span("Future work")),
        tags$div(class = "module-state pending", tags$strong("Interactive dashboard"), tags$span("Scaffold only"))
      )
    )
  )
)

server <- function(input, output, session) {
  output$summary_cards <- renderUI({
    figure_outputs_dir <- file.path(repo_root, "manuscript", "figures", "outputs")
    main_dir <- file.path(figure_outputs_dir, "main")
    supplementary_figures_dir <- file.path(figure_outputs_dir, "extended_data")
    legacy_dir <- file.path(figure_outputs_dir, "legacy")
    figure_data_dir <- file.path(repo_root, "manuscript", "figure_data")

    cards <- list(
      list(label = "Main figures", value = count_assets(main_dir), tone = "warm"),
      list(label = "Supplementary figures", value = count_assets(supplementary_figures_dir), tone = "cool"),
      list(label = "Legacy exports", value = count_assets(legacy_dir), tone = "muted"),
      list(label = "Figure-data extracts", value = count_assets(figure_data_dir, pattern = "\\.tsv$"), tone = "forest")
    )

    div(
      class = "summary-grid",
      lapply(cards, function(card) {
        div(
          class = paste("summary-card", card$tone),
          tags$div(class = "summary-value", card$value),
          tags$div(class = "summary-label", card$label)
        )
      })
    )
  })

  output$asset_table <- renderTable({
    figure_outputs_dir <- file.path(repo_root, "manuscript", "figures", "outputs")
    rbind(
      list_assets(file.path(figure_outputs_dir, "main"), "main"),
      list_assets(file.path(figure_outputs_dir, "extended_data"), "supplementary_figures"),
      list_assets(file.path(figure_outputs_dir, "legacy"), "legacy")
    )
  }, striped = TRUE, bordered = FALSE, spacing = "s")

  output$figure_data_table <- renderTable({
    list_assets(file.path(repo_root, "manuscript", "figure_data"), "figure_data")
  }, striped = TRUE, bordered = FALSE, spacing = "s")
}

shinyApp(ui = ui, server = server)
