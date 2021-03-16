####################################################
####################################################
drsimonj_colors <- c(`red` = "#d11141", `green` = "#00b159", `blue` = "#00aedb", `orange` = "#f37735", 
                     `yellow` = "#ffc425", `light grey` = "#cccccc", `dark grey` = "#8c8c8c")
####################################################
####################################################
drsimonj_cols <- function(...) {
   cols <- c(...)
   if (is.null(cols))
     return (drsimonj_colors)
   drsimonj_colors[cols]
 }

####################################################
####################################################
custom_palettes <- c(wesanderson::wes_palettes, 
                     list(`main`  = drsimonj_cols("blue", "green", "yellow"),
  `cool`  = drsimonj_cols("blue", "green"),
  `hot`   = drsimonj_cols("yellow", "orange", "red"),
  `mixed` = drsimonj_cols("blue", "green", "yellow", "orange", "red"),
  `grey`  = drsimonj_cols("light grey", "dark grey"),
  `cbp1` = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
  `cbp2` = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
  `Partners1` = c("#FF0F1D", "#0095C8", "#00BB56", "#BD64B0", "#FF8A00", "#D16626"),
  `Partners2` = c("#FF2E2B", "#0095CF", "#00C156", "#FFAD00", "#CB80BB", "#CBC9C9"),
  `Set1` = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999"),
  `Set2` = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
  `Set3` = c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
  `Dark2` = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
  `Accent` = c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666"),
  `Paired` = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
  `Pastel1` = c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2"),
  `Pastel2` = c("#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC"),
  `Projects` = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
  `Clusters` = c("#bd64b0", "#ff8a00", "#00bb56", "#0095c8", "#ff0f1d"))
)

####################################################
####################################################

#' Return function to interpolate a drsimonj color palette
#'
#' @param palette Character name of palette in custom_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette()
#'
drsimonj_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- custom_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}


####################################################
####################################################

#' Fill scale constructor for drsimonj colors
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_drsimonj <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- drsimonj_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("drsimonj_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}

