suppressPackageStartupMessages({ library(ggplot2) })
plot_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 12)
    )
}
write_plot <- function(p, stem, width = 3000, height = 2200, res = 300, pdf = FALSE,
                       width_in = 10, height_in = 8) {
  tiff(file.path("figures", paste0(stem, ".tiff")), width = width, height = height, res = res)
  print(p); dev.off()
  if (pdf) {
    pdf(file.path("figures", paste0(stem, ".pdf")), width = width_in, height = height_in)
    print(p); dev.off()
  }
}
