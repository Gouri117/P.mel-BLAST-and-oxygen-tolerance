library(ComplexHeatmap)
library(circlize)
library(grid)

matrix_file <- "species_gene_ppos_matrix.csv"
out_pdf     <- "bacteria_dendrogram_heatmap_compact.pdf"

## Read matrix
df <- read.csv(matrix_file, row.names = 1, check.names = FALSE)
mat <- as.matrix(df)
mat <- apply(mat, 2, function(x) as.numeric(as.character(x)))
rownames(mat) <- rownames(df)
mat[is.na(mat)] <- 0

## Cluster bacteria (rows)
hc <- hclust(dist(mat), method = "complete")

## Transpose: rows = proteins, cols = bacteria
mat_plot <- t(mat)

## Custom red scale (as requested)
col_fun <- colorRamp2(
  c(0, 40, 60, 80, 100),
  c(
    "lightyellow",
    "orange1",   # very light red
    "orange2",   # light red
    "orange3",   # red
    "#cb181d"    # dark red
  )
)

## Compact PDF
pdf(out_pdf, width = 14, height = 8)

ht <- Heatmap(
  mat_plot,
  name = "PPOS",
  col = col_fun,
  
  ## Bacteria dendrogram
  cluster_columns = as.dendrogram(hc),
  show_column_dend = TRUE,
  column_dend_height = unit(3.5, "cm"),
  
  ## Make boxes smaller (key change)
  width  = unit(ncol(mat_plot) * 2.5, "mm"),
  height = unit(nrow(mat_plot) * 2.5, "mm"),
  
  ## Proteins
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  
  ## Species names: italic + compact
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_gp = gpar(
    fontsize = 8,
    fontface = "italic"
  ),
  
  heatmap_legend_param = list(
    title = "PPOS",
    at = c(0, 40, 60, 80, 100)
  )
)

draw(ht, heatmap_legend_side = "right")
dev.off()

message("Wrote: ", out_pdf)


















