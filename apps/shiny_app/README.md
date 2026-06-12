# Interactive Visualization Scaffold

This directory contains the repository's interactive delivery scaffold.

## Current Status

- The app is a lightweight Shiny workspace, not a production dashboard.
- It surfaces the current figure inventory and manuscript-facing assets.
- It is intended to sit on top of curated outputs, not raw pipeline intermediates.

## Run

```bash
R -e "shiny::runApp('shiny_app')"
```

## Design Goals

- Reuse the same visual language as the static figure workspace.
- Make the current asset inventory visible to collaborators.
- Keep future interactive modules separate from the static publication pipeline.
