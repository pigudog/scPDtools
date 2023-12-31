% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCP-plot.R
\name{palette_scp}
\alias{palette_scp}
\title{Color palettes collected in SCP.}
\usage{
palette_scp(
  x,
  n = 100,
  palette = "Paired",
  palcolor = NULL,
  type = "auto",
  matched = FALSE,
  reverse = FALSE,
  NA_keep = FALSE,
  NA_color = "grey80"
)
}
\arguments{
\item{x}{A vector of character/factor or numeric values. If missing, numeric values 1:n will be used as x.}

\item{n}{The number of colors to return for numeric values.}

\item{palette}{Palette name. All available palette names can be queried with \code{show_palettes()}.}

\item{palcolor}{Custom colors used to create a color palette.}

\item{type}{Type of \code{x}. Can be one of "auto", "discrete" or "continuous". The default is "auto", which automatically detects if \code{x} is a numeric value.}

\item{matched}{If \code{TRUE}, will return a color vector of the same length as \code{x}.}

\item{reverse}{Whether to invert the colors.}

\item{NA_keep}{Whether to keep the color assignment to NA in \code{x}.}

\item{NA_color}{Color assigned to NA if NA_keep is \code{TRUE}.}
}
\description{
Color palettes collected in SCP.
}
\examples{
x <- c(1:3, NA, 3:5)
(pal1 <- palette_scp(x, palette = "Spectral"))
(pal2 <- palette_scp(x, palcolor = c("red", "white", "blue")))
(pal3 <- palette_scp(x, palette = "Spectral", n = 10))
(pal4 <- palette_scp(x, palette = "Spectral", n = 10, reverse = TRUE))
(pal5 <- palette_scp(x, palette = "Spectral", matched = TRUE))
(pal6 <- palette_scp(x, palette = "Spectral", matched = TRUE, NA_keep = TRUE))
(pal7 <- palette_scp(x, palette = "Paired", type = "discrete"))
show_palettes(list(pal1, pal2, pal3, pal4, pal5, pal6, pal7))

all_palettes <- show_palettes(return_palettes = TRUE)
names(all_palettes)

\dontrun{
if (interactive()) {
  check_R(c("stringr", "RColorBrewer", "ggsci", "Redmonder", "rcartocolor", "nord", "viridis", "pals", "oompaBase", "dichromat", "jaredhuling/jcolors"))
  library(stringr)
  library(RColorBrewer)
  library(ggsci)
  library(Redmonder)
  library(rcartocolor)
  library(nord)
  library(viridis)
  library(pals)
  library(dichromat)
  library(jcolors)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  ggsci_db <- ggsci:::ggsci_db
  redmonder.pal.info <- Redmonder::redmonder.pal.info
  metacartocolors <- rcartocolor::metacartocolors
  rownames(metacartocolors) <- metacartocolors$Name
  nord_palettes <- nord::nord_palettes
  viridis_names <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
  viridis_palettes <- lapply(setNames(viridis_names, viridis_names), function(x) viridis::viridis(100, option = x))
  ocean_names <- names(pals:::syspals)[grep("ocean", names(pals:::syspals))]
  ocean_palettes <- pals:::syspals[ocean_names]
  dichromat_palettes <- dichromat::colorschemes
  jcolors_names <- paste0("jcolors-", c("default", "pal2", "pal3", "pal4", "pal5", "pal6", "pal7", "pal8", "pal9", "pal10", "pal11", "pal12", "rainbow"))
  custom_names <- c("jet", "simspec", "GdRd")
  custom_palettes <- list(
    oompaBase::jetColors(N = 100),
    c("#c22b86", "#f769a1", "#fcc5c1", "#253777", "#1d92c0", "#9ec9e1", "#015b33", "#42aa5e", "#d9f0a2", "#E66F00", "#f18c28", "#FFBB61"),
    c("gold", "red3")
  )
  names(custom_palettes) <- custom_names

  palette_list <- list()
  all_colors <- c(
    rownames(brewer.pal.info), names(ggsci_db), rownames(redmonder.pal.info),
    rownames(metacartocolors), names(nord_palettes), names(viridis_palettes),
    ocean_names, names(dichromat_palettes), jcolors_names,
    custom_names
  )
  for (pal in all_colors) {
    if (!pal \%in\% all_colors) {
      stop(paste0("Invalid pal Must be one of ", paste0(all_colors, collapse = ",")))
    }
    if (pal \%in\% rownames(brewer.pal.info)) {
      pal_n <- brewer.pal.info[pal, "maxcolors"]
      pal_category <- brewer.pal.info[pal, "category"]
      if (pal_category == "div") {
        palcolor <- rev(brewer.pal(name = pal, n = pal_n))
      } else {
        if (pal == "Paired") {
          palcolor <- brewer.pal(12, "Paired")[c(1:4, 7, 8, 5, 6, 9, 10, 11, 12)]
        } else {
          palcolor <- brewer.pal(name = pal, n = pal_n)
        }
      }
      if (pal_category == "qual") {
        attr(palcolor, "type") <- "discrete"
      } else {
        attr(palcolor, "type") <- "continuous"
      }
    } else if (pal \%in\% names(ggsci_db)) {
      if (pal \%in\% c("d3", "uchicago", "material")) {
        for (subpal in names(ggsci_db[[pal]])) {
          palcolor <- ggsci_db[[pal]][[subpal]]
          if (pal == "material") {
            attr(palcolor, "type") <- "continuous"
          } else {
            attr(palcolor, "type") <- "discrete"
          }
          palette_list[[paste0(pal, "-", subpal)]] <- palcolor
        }
        next
      } else {
        palcolor <- ggsci_db[[pal]][[1]]
        if (pal == "gsea") {
          attr(palcolor, "type") <- "continuous"
        } else {
          attr(palcolor, "type") <- "discrete"
        }
      }
    } else if (pal \%in\% rownames(redmonder.pal.info)) {
      pal_n <- redmonder.pal.info[pal, "maxcolors"]
      pal_category <- redmonder.pal.info[pal, "category"]
      if (pal_category == "div") {
        palcolor <- rev(redmonder.pal(name = pal, n = pal_n))
      } else {
        palcolor <- redmonder.pal(name = pal, n = pal_n)
      }
      if (pal_category == "qual") {
        attr(palcolor, "type") <- "discrete"
      } else {
        attr(palcolor, "type") <- "continuous"
      }
    } else if (pal \%in\% rownames(metacartocolors)) {
      pal_n <- metacartocolors[pal, "Max_n"]
      palcolor <- carto_pal(name = pal, n = pal_n)
      if (pal_category == "qualitative") {
        attr(palcolor, "type") <- "discrete"
      } else {
        attr(palcolor, "type") <- "continuous"
      }
    } else if (pal \%in\% names(nord_palettes)) {
      palcolor <- nord_palettes[[pal]]
      attr(palcolor, "type") <- "discrete"
    } else if (pal \%in\% names(viridis_palettes)) {
      palcolor <- viridis_palettes[[pal]]
      attr(palcolor, "type") <- "continuous"
    } else if (pal \%in\% names(ocean_palettes)) {
      palcolor <- ocean_palettes[[pal]]
      attr(palcolor, "type") <- "continuous"
    } else if (pal \%in\% names(dichromat_palettes)) {
      palcolor <- dichromat_palettes[[pal]]
      if (pal \%in\% c("Categorical.12", "SteppedSequential.5")) {
        attr(palcolor, "type") <- "discrete"
      } else {
        attr(palcolor, "type") <- "continuous"
      }
    } else if (pal \%in\% jcolors_names) {
      palcolor <- jcolors(palette = gsub("jcolors-", "", pal))
      if (pal \%in\% paste0("jcolors-", c("pal10", "pal11", "pal12", "rainbow"))) {
        attr(palcolor, "type") <- "continuous"
      } else {
        attr(palcolor, "type") <- "discrete"
      }
    } else if (pal \%in\% custom_names) {
      palcolor <- custom_palettes[[pal]]
      if (pal \%in\% c("jet")) {
        attr(palcolor, "type") <- "continuous"
      } else {
        attr(palcolor, "type") <- "discrete"
      }
    }
    palette_list[[pal]] <- palcolor
  }
  # usethis::use_data(palette_list, internal = TRUE)
}
}
}
\seealso{
\code{\link{show_palettes}}
}
