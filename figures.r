suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(infercnv)
  library(forcats)
})

SO.SEQ <- readRDS("out/seurat/dropseq.seurat.rds")
SO.SPATIAL <- readRDS("out/seurat/merscope.seurat.rds")
TMA.COORDS <- fread("data/merscope/transcripts/tma-origin.csv")
LM548.COORDS <- fread("data/merscope/transcripts/lm548-origin.csv")
PANEL.META <- fread("data/merscope/panel-meta.csv")[, .(
  gene = `Vizgen Gene`,
  module = str_replace_all(Notes, " |:|/", "_") |>
    str_remove_all("\\(|\\)")
)]

DURA.LM548.COORDS.X <- c(3700, 4700)
DURA.LM548.COORDS.Y <- c(1400, 2100)
DURA.LM548.COORDS.THETA <- -30
BRAIN.LM548.COORDS.X <- c(3200, 4200)
BRAIN.LM548.COORDS.Y <- c(-2725, -2425)
BRAIN.LM548.COORDS.THETA <- -110
DURA.LM539.COORDS.X <- c(-500, 500)
DURA.LM539.COORDS.Y <- c(6800, 7500)
DURA.LM539.COORDS.THETA <- 20
BRAIN.LM539.COORDS.X <- c(-12200, -11200)
BRAIN.LM539.COORDS.Y <- c(5525, 5825)
BRAIN.LM539.COORDS.THETA <- 100

PLOT.CELL.POINT.SIZE <- 1.5
PLOT.SCALE.BAR.SIZE <- 100

PAL.FIG2B <- c(
  "dura" = "#3283fe",
  "arachnoid" = "#0ba71d",
  "pia" = "#fae107",
  "mural" = "#e99a12",
  "endothelium" = "#f6222e",
  "myeloid" = "#fe00fa",
  "leukocyte" = "#e0a8cd"
)
PAL.FIG2C <- c(
  "dura" = "#3283fe",
  "arachnoid" = "#0ba71d",
  "pia" = "#fae107",
  "mural" = "#e99a12",
  "endothelium" = "#f6222e",
  "myeloid" = "#fe00fa",
  "leukocyte" = "#e0a8cd",
  "parenchyma" = "#f0ead6",
  "unknown" = "#B00068"
)
PAL.FIG6E <- c(
  "IGF2" = "#6a1c89",
  "IGF1R" = "#f07819",
  "IGFBP7" = "#07bdf3",
  "dura" = "#b3b5b8",
  "arachnoid" = "#b3b5b8",
  "pia" = "#b3b5b8",
  "mural" = "#b3b5b8",
  "endothelium" = "#b3b5b8"
)

cdt <- function(DT, theta, dtm) {
  copy(DT)[, `:=`(
    rot.x = global_x * cos(theta) - global_y * sin(theta),
    rot.y = global_x * sin(theta) + global_y * cos(theta)
  )][
    dtm,
    on = .(gene)
  ][, .(rot.x, rot.y, gene, module)]
}

gdc <- function(so, col) {
  if (col %chin% colnames(so[[]])) {
    as.data.table(so[[col]], keep.rownames = "cell")
  } else {
    GetAssayData(so)[col, ] |>
      as.data.table(keep.rownames = TRUE) |>
      `colnames<-`(c("cell", col))
  }
}

sdp <- function(so, fov, theta, x.lims, y.lims, size, fill.col, outline.col = "#00000000", stroke.size = 0, pal.v, lims.v, pal.m) {
  po <- as.data.table(GetTissueCoordinates(so, image = fov))[
    gdc(so, fill.col),
    on = .(cell)
  ][, .(
    rot.x = x * cos(theta) - y * sin(theta),
    rot.y = x * sin(theta) + y * cos(theta),
    fill.col = get(fill.col)
  )][
    rot.x %between% x.lims & rot.y %between% y.lims
  ] |>
    ggplot() +
    aes(
      x = rot.x,
      y = rot.y
    ) +
    geom_point(size = size, shape = 21, stroke = stroke.size, colour = outline.col) +
    coord_fixed()

  if (!missing(pal.m)) {
    po <- po + scale_fill_manual(values = pal.m) + aes(fill = factor(fill.col, levels = names(pal.m)))
  } else {
    po <- po + aes(fill = fill.col)
  }
  if (!missing(pal.v)) {
    po <- po + scale_fill_viridis_c(option = pal.v, limits = lims.v)
  }
  po
}

tdp <- function(coords, include.genes, include.modules, x.lims, y.lims, alpha.def = 1, palette = NULL, pt.size = NULL) {
  conv <- function(l) {
    list(
      if (is.null(names(l))) {
        l
      } else {
        names(l)
      },
      l
    ) |>
      pmap(function(n, v) {
        if (is.numeric(v)) {
          data.table(alpha = v, col = n)
        } else {
          data.table(alpha = alpha.def, col = v)
        }
      }) |>
      reduce(rbind, .init = data.table(col = character(), alpha = numeric()))
  }

  include.genes <- conv(include.genes)
  include.modules <- conv(include.modules)

  po <- rbind(
    coords[
      rot.x %between% x.lims & rot.y %between% y.lims,
      .(
        col = gene,
        rot.x,
        rot.y
      )
    ][include.genes, on = .(col)],
    coords[
      rot.x %between% x.lims & rot.y %between% y.lims,
      .(
        col = module,
        rot.x,
        rot.y
      )
    ][include.modules, on = .(col)]
  ) |>
    ggplot() +
    aes(x = rot.x, y = rot.y, alpha = alpha) +
    coord_fixed() +
    scale_alpha_identity() +
    guides(
      alpha = "none"
    )

  if (is.null(pt.size)) {
    po <- po + geom_point(shape = ".")
  } else {
    po <- po + geom_point(shape = 21, size = pt.size)
  }

  if (is.null(palette)) {
    po <- po +
      scale_colour_manual(
        name = NULL,
        values = ebp("Set1", length(c(include.genes, include.modules)))
      ) +
      aes(colour = factor(col))
  } else {
    po <- po +
      scale_colour_manual(
        name = NULL,
        values = palette
      ) +
      aes(colour = factor(col, levels = names(palette)))
  }
  po
}

ebp <- function(palette, n) {
  if (n > brewer.pal.info[palette, "maxcolors"]) {
    colorRampPalette(
      brewer.pal(
        brewer.pal.info[palette, "maxcolors"],
        palette
      )
    )(n)
  } else {
    brewer.pal(
      brewer.pal.info[palette, "maxcolors"],
      palette
    )
  }
}

# figure 2b
o <- wrap_plots(
  tdp(
    cdt(TMA.COORDS, DURA.LM548.COORDS.THETA * pi / 180, PANEL.META)[!gene %chin% c("ATP5F1E", "ATP5MG", "ATP1A2")],
    c(), c("dura", "arachnoid", "pia", "endothelium", "mural", "leukocyte", "myeloid"),
    DURA.LM548.COORDS.X, DURA.LM548.COORDS.Y,
    0.5,
    PAL.FIG2B
  ),
  tdp(
    cdt(LM548.COORDS, BRAIN.LM548.COORDS.THETA * pi / 180, PANEL.META)[!gene %chin% c("ATP5F1E", "ATP5MG", "ATP1A2")],
    c(), c("dura", "arachnoid", "pia", "endothelium", "mural", "leukocyte", "myeloid"),
    BRAIN.LM548.COORDS.X, BRAIN.LM548.COORDS.Y,
    0.5,
    PAL.FIG2B
  ) + geom_path(aes(x = x, y = y), data.frame(x = c(BRAIN.LM548.COORDS.X[[2]] - PLOT.SCALE.BAR.SIZE, BRAIN.LM548.COORDS.X[[2]]), y = rep(BRAIN.LM548.COORDS.Y[[1]] - 50, 2)), colour = "red", inherit.aes = FALSE),
  nrow = 2,
  guides = "collect"
) &
  theme_void() &
  guides(
    color = guide_legend(
      override.aes = list(size = 3, shape = "circle"),
      title = NULL
    )
  ) &
  theme(
    text = element_text(size = 26, colour = "white"),
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave(
  "out/fig-2b-lm548.png", o,
  w = 4000, h = 4000, units = "px", create.dir = TRUE
)
o <- wrap_plots(
  tdp(
    cdt(TMA.COORDS, DURA.LM539.COORDS.THETA * pi / 180, PANEL.META)[!gene %chin% c("ATP5F1E", "ATP5MG", "ATP1A2")],
    c(), c("pia", "dura", "arachnoid", "endothelium", "mural", "leukocyte", "myeloid"),
    DURA.LM539.COORDS.X, DURA.LM539.COORDS.Y,
    0.5,
    PAL.FIG2B
  ) + scale_y_reverse(),
  tdp(
    cdt(TMA.COORDS, BRAIN.LM539.COORDS.THETA * pi / 180, PANEL.META)[!gene %chin% c("ATP5F1E", "ATP5MG", "ATP1A2")],
    c(), c("pia", "dura", "arachnoid", "endothelium", "mural", "leukocyte", "myeloid"),
    BRAIN.LM539.COORDS.X, BRAIN.LM539.COORDS.Y,
    0.5,
    PAL.FIG2B
  ) + geom_path(aes(x = x, y = y), data.frame(x = c(BRAIN.LM539.COORDS.X[[2]] - PLOT.SCALE.BAR.SIZE, BRAIN.LM539.COORDS.X[[2]]), y = rep(BRAIN.LM539.COORDS.Y[[1]] - 50, 2)), colour = "red", inherit.aes = FALSE),
  nrow = 2,
  guides = "collect"
) &
  theme_void() &
  guides(
    color = guide_legend(
      override.aes = list(size = 3, shape = "circle"),
      title = NULL
    )
  ) &
  theme(
    text = element_text(size = 26, colour = "white"),
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave(
  "out/fig-2b-lm539.png", o,
  w = 4000, h = 4000, units = "px", create.dir = TRUE
)

# figure 2c
o <- wrap_plots(
  sdp(
    SO.SPATIAL,
    "tma.cnt",
    DURA.LM548.COORDS.THETA * pi / 180,
    DURA.LM548.COORDS.X, DURA.LM548.COORDS.Y,
    PLOT.CELL.POINT.SIZE,
    "cell_type",
    pal.m = PAL.FIG2C
  ) +
    labs(fill = NULL, x = NULL, y = NULL),
  sdp(
    SO.SPATIAL,
    "lm548.cnt",
    BRAIN.LM548.COORDS.THETA * pi / 180,
    BRAIN.LM548.COORDS.X, BRAIN.LM548.COORDS.Y,
    PLOT.CELL.POINT.SIZE,
    "cell_type",
    pal.m = PAL.FIG2C
  ) +
    labs(fill = NULL, x = NULL, y = NULL) +
    geom_path(aes(x = x, y = y), data.frame(x = c(BRAIN.LM548.COORDS.X[[2]] - PLOT.SCALE.BAR.SIZE, BRAIN.LM548.COORDS.X[[2]]), y = rep(BRAIN.LM548.COORDS.Y[[1]] - 50, 2)), colour = "red", inherit.aes = FALSE),
  guides = "collect",
  nrow = 2
) &
  guides(
    color = guide_legend(
      override.aes = list(size = 3, shape = "circle"),
      title = NULL
    )
  ) &
  theme_void() &
  theme(
    axis.text = element_text(size = 0.25),
    axis.ticks = element_line(),
    axis.title = element_blank(),
    text = element_text(size = 26, colour = "white"),
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave(
  "out/fig-2c-lm548.png",
  plot = o,
  w = 3000, h = 3000, units = "px", create.dir = TRUE
)

o <- wrap_plots(
  sdp(
    SO.SPATIAL,
    "tma.cnt",
    DURA.LM539.COORDS.THETA * pi / 180,
    DURA.LM539.COORDS.X, DURA.LM539.COORDS.Y,
    PLOT.CELL.POINT.SIZE,
    "cell_type",
    pal.m = PAL.FIG2C
  ) +
    labs(
      fill = NULL, x = NULL, y = NULL
    ) + scale_y_reverse(),
  sdp(
    SO.SPATIAL,
    "tma.cnt",
    BRAIN.LM539.COORDS.THETA * pi / 180,
    BRAIN.LM539.COORDS.X, BRAIN.LM539.COORDS.Y,
    PLOT.CELL.POINT.SIZE,
    "cell_type",
    pal.m = PAL.FIG2C
  ) +
    labs(fill = NULL, x = NULL, y = NULL) +
    geom_path(aes(x = x, y = y), data.frame(x = c(BRAIN.LM539.COORDS.X[[2]] - PLOT.SCALE.BAR.SIZE, BRAIN.LM539.COORDS.X[[2]]), y = rep(BRAIN.LM539.COORDS.Y[[1]] - 50, 2)), colour = "red", inherit.aes = FALSE),
  guides = "collect",
  nrow = 2
) &
  guides(
    color = guide_legend(
      override.aes = list(size = 3, shape = "circle"),
      title = NULL
    )
  ) &
  theme_void() &
  theme(
    axis.text = element_text(size = 0.25),
    axis.ticks = element_line(),
    axis.title = element_blank(),
    text = element_text(size = 26, colour = "white"),
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave(
  "out/fig-2c-lm539.png",
  plot = o,
  w = 3000, h = 3000, units = "px", create.dir = TRUE
)

# figure 6e
o <- wrap_plots(
  tdp(
    cdt(TMA.COORDS, DURA.LM548.COORDS.THETA * pi / 180, PANEL.META)[!gene %chin% c("ATP5F1E", "ATP5MG", "ATP1A2")],
    c("IGF2" = 1, "IGF1R" = 1, "IGFBP7" = 1), c("endothelium", "arachnoid", "dura", "pia", "mural"),
    DURA.LM548.COORDS.X, DURA.LM548.COORDS.Y,
    0.3,
    PAL.FIG6E
  ),
  tdp(
    cdt(LM548.COORDS, BRAIN.LM548.COORDS.THETA * pi / 180, PANEL.META)[!gene %chin% c("ATP5F1E", "ATP5MG", "ATP1A2")],
    c("IGF2" = 1, "IGF1R" = 1, "IGFBP7" = 1), c("endothelium", "arachnoid", "dura", "pia", "mural"),
    BRAIN.LM548.COORDS.X, BRAIN.LM548.COORDS.Y,
    0.3,
    PAL.FIG6E
  ) + geom_path(aes(x = x, y = y), data.frame(x = c(BRAIN.LM548.COORDS.X[[2]] - PLOT.SCALE.BAR.SIZE, BRAIN.LM548.COORDS.X[[2]]), y = rep(BRAIN.LM548.COORDS.Y[[1]] - 50, 2)), colour = "red", inherit.aes = FALSE),
  nrow = 2,
  guides = "collect"
) &
  theme_void() &
  guides(
    color = guide_legend(
      override.aes = list(size = 3, shape = "circle"),
      title = NULL
    )
  ) &
  theme(
    text = element_text(size = 26),
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave(
  "out/fig-6e-lm548.png", o,
  w = 4000, h = 4000, units = "px", create.dir = TRUE
)
o <- wrap_plots(
  tdp(
    cdt(TMA.COORDS, DURA.LM539.COORDS.THETA * pi / 180, PANEL.META)[!gene %chin% c("ATP5F1E", "ATP5MG", "ATP1A2")],
    c("IGF2" = 1, "IGF1R" = 1, "IGFBP7" = 1), c("endothelium", "arachnoid", "dura", "pia", "mural"),
    DURA.LM539.COORDS.X, DURA.LM539.COORDS.Y,
    0.3,
    PAL.FIG6E
  ) + scale_y_reverse(),
  tdp(
    cdt(TMA.COORDS, BRAIN.LM539.COORDS.THETA * pi / 180, PANEL.META)[!gene %chin% c("ATP5F1E", "ATP5MG", "ATP1A2")],
    c("IGF2" = 1, "IGF1R" = 1, "IGFBP7" = 1), c("endothelium", "arachnoid", "dura", "pia", "mural"),
    BRAIN.LM539.COORDS.X, BRAIN.LM539.COORDS.Y,
    0.3,
    PAL.FIG6E
  ) + geom_path(aes(x = x, y = y), data.frame(x = c(BRAIN.LM539.COORDS.X[[2]] - PLOT.SCALE.BAR.SIZE, BRAIN.LM539.COORDS.X[[2]]), y = rep(BRAIN.LM539.COORDS.Y[[1]] - 50, 2)), colour = "red", inherit.aes = FALSE),
  nrow = 2,
  guides = "collect"
) &
  theme_void() &
  guides(
    color = guide_legend(
      override.aes = list(size = 3, shape = "circle"),
      title = NULL
    )
  ) &
  theme(
    text = element_text(size = 26),
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave(
  "out/fig-6e-lm539.png", o,
  w = 4000, h = 4000, units = "px", create.dir = TRUE
)

# supplementary figure 1b
ICNV.META <- as.data.table(SO.SEQ[[c("cell_type", "patient")]], keep.rownames = "cell")[, .(
  cell,
  patient,
  icnv.col = pmap_vec(list(cell_type, patient), function(cell_type, patient) {
    if (str_detect(cell_type, "Myeloid|T/NK|Neutrophils|Macrophages|Bcells|Mast")) {
      "ref"
    } else if (str_detect(cell_type, "AraC[0-9]")) {
      sprintf("SA_%s", patient)
    } else if (str_detect(cell_type, "PiaC[0-9]")) {
      sprintf("SP_%s", patient)
    } else if (str_detect(cell_type, "DuraC[0-9]")) {
      sprintf("SD_%s", patient)
    } else {
      NA
    }
  })
)][!is.na(icnv.col)]
ICNV.META <- rbind(
  # mix together 100 immune cells from each patient
  ICNV.META[icnv.col == "ref", first(.SD[sample(.N)], 100), by = .(patient)],
  # mix together 300 fibroblast cells from each patient
  ICNV.META[icnv.col != "ref", first(.SD[sample(.N)], 300), by = .(patient)]
)
CreateInfercnvObject(
  GetAssayData(
    subset(SO.SEQ, cells = ICNV.META[, cell]),
    "RNA", "counts"
  ),
  data.frame(
    fread("https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt", sep = "\t", fill = TRUE)[, 1:4],
    row.names = "V1"
  ),
  data.frame(ICNV.META[, !c("patient")], row.names = "cell"),
  "ref"
) |>
  run(
    cutoff = 0.1,
    out_dir = tempdir(),
    cluster_by_groups = TRUE, denoise = TRUE, HMM = TRUE,
    num_threads = 1, resume_mode = FALSE, no_plot = TRUE
  ) |>
  plot_cnv(contig_cex = 2.1, out_dir = "out/figs1b")

# supplementary figure 5
walk(
  c("ALCAM", "APOD", "MGP"),
  function(col_oi) {
    lims <- c(
      gdc(SO.SPATIAL, col_oi)[, min(get(col_oi))],
      gdc(SO.SPATIAL, col_oi)[, max(get(col_oi))]
    )
    o <- wrap_plots(
      sdp(
        SO.SPATIAL,
        "tma.cnt",
        DURA.LM548.COORDS.THETA * pi / 180,
        DURA.LM548.COORDS.X, DURA.LM548.COORDS.Y,
        PLOT.CELL.POINT.SIZE,
        col_oi,
        pal.v = "turbo",
        lims.v = lims
      ),
      sdp(
        SO.SPATIAL,
        "lm548.cnt",
        BRAIN.LM548.COORDS.THETA * pi / 180,
        BRAIN.LM548.COORDS.X, BRAIN.LM548.COORDS.Y,
        PLOT.CELL.POINT.SIZE,
        col_oi,
        pal.v = "turbo",
        lims.v = lims
      ) +
        geom_path(aes(x = x, y = y), data.frame(x = c(BRAIN.LM548.COORDS.X[[2]] - PLOT.SCALE.BAR.SIZE, BRAIN.LM548.COORDS.X[[2]]), y = rep(BRAIN.LM548.COORDS.Y[[1]] - 50, 2)), colour = "red", inherit.aes = FALSE),
      guides = "collect",
      nrow = 2
    ) &
      theme_void() &
      theme(
        axis.text = element_text(size = 0.25),
        axis.ticks = element_line(),
        axis.title = element_blank(),
        text = element_text(size = 26, colour = "white"),
        plot.margin = margin(0, 0, 0, 0)
      ) &
      labs(colour = NULL, x = NULL, y = NULL, fill = col_oi)
    ggsave(
      sprintf("out/fig-s6-lm548-%s.png", col_oi),
      plot = o,
      w = 2500, h = 2500, units = "px", create.dir = TRUE
    )
    o <- wrap_plots(
      sdp(
        SO.SPATIAL,
        "tma.cnt",
        DURA.LM539.COORDS.THETA * pi / 180,
        DURA.LM539.COORDS.X, DURA.LM539.COORDS.Y,
        PLOT.CELL.POINT.SIZE,
        col_oi,
        pal.v = "turbo",
        lims.v = lims
      ) + scale_y_reverse(),
      sdp(
        SO.SPATIAL,
        "tma.cnt",
        BRAIN.LM539.COORDS.THETA * pi / 180,
        BRAIN.LM539.COORDS.X, BRAIN.LM539.COORDS.Y,
        PLOT.CELL.POINT.SIZE,
        col_oi,
        pal.v = "turbo",
        lims.v = lims
      ) +
        geom_path(aes(x = x, y = y), data.frame(x = c(BRAIN.LM539.COORDS.X[[2]] - PLOT.SCALE.BAR.SIZE, BRAIN.LM539.COORDS.X[[2]]), y = rep(BRAIN.LM539.COORDS.Y[[1]] - 50, 2)), colour = "red", inherit.aes = FALSE),
      guides = "collect",
      nrow = 2
    ) &
      theme_void() &
      theme(
        axis.text = element_text(size = 0.25),
        axis.ticks = element_line(),
        axis.title = element_blank(),
        text = element_text(size = 26, colour = "white"),
        plot.margin = margin(0, 0, 0, 0)
      ) &
      labs(colour = NULL, x = NULL, y = NULL, fill = col_oi)
    ggsave(
      sprintf("out/fig-s6-lm539-%s.png", col_oi),
      plot = o,
      w = 2500, h = 2500, units = "px", create.dir = TRUE
    )
  }
)
