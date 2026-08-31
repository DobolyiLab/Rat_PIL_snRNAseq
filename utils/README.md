# utils/

Place your `plotting.R` here. It is not included in this repository -- you
already have it (see the DobolyiLab conversation this repo came out of) --
but every script that produces figures expects it to define two functions
with this interface:

```r
plot_theme()
# Returns a ggplot2 theme (or list of theme elements) to be added to a plot
# with `+`, e.g. DimPlot(obj) + plot_theme().

write_plot(plot, name, pdf = FALSE)
# Writes `plot` to figures/<name>.tiff (matching journal TIFF requirements),
# and additionally to figures/<name>.pdf when pdf = TRUE.
```

`scripts/celltype_labeling.R` and `scripts/Fos_Calb1_quantification.R` call
`source("utils/plotting.R")` and use both functions with this signature. If
your actual file uses different function names or argument order, either
rename to match, or edit the two calling scripts to match your file --
whichever is less disruptive to the rest of your codebase.
