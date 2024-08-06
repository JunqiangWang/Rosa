#' @description This function does the waterfall plot for a set of genes. The x-axis is the ranked signature and the y-axis
#' is the Z scores transformed by the ranked signature
#'
#' @export
WaterfallGSEA2 <- function(
    signature = signature,
    geneset = geneset,
    gsea.out = gsea.out,
    show.labels = TRUE,
    label.cutoff.z = 0.5
) {

  require(graphics)

  # Ensure Arial font is available
 # if (!("Arial" %in% fonts())) {
   # warning("Arial font not found. Using default font.")
  #}

  # Set font family to Arial

  signature <- sort(signature, decreasing = TRUE)
  rank.signature <- rank(signature)

  ledge <- intersect(gsea.out$ledge %>% names(), names(signature[signature > label.cutoff.z]))

  df <- data.frame(gene = names(rank.signature), Stouffer.signature = signature, rank.signature = rank.signature)
  df$label <- ifelse(df$gene %in% ledge, 1, 0)
  df$label2 <- ifelse(df$label == 1, df$gene, NA)

  # Create a basic scatter plot using base R plot function
  # Create a basic scatter plot using base R plot function
  plot(df$rank.signature, df$Stouffer.signature, col = "grey", pch = 16,
       xlim = c(min(df$rank.signature), max(df$rank.signature) * 2.2),
       xlab = "Rank Signature", ylab = "Differential Signature",
       main = "Waterfall Plot")

  # Overlay the red dots on top
  points(df$rank.signature[df$label == 1], df$Stouffer.signature[df$label == 1], col = "orange", pch = 16)


  if (show.labels && length(which(df$label == 1)) > 0) {
    s11 <- which(df$label == 1)

    # Calculate vertical spacing for labels
    dd <- (par("usr")[4] - par("usr")[3]) / (length(s11) + 1)
    ypos <- seq(from = par("usr")[3], by = dd, length.out = length(s11)) + dd
    ypos <- rev(ypos)

    # Add labels and connecting segments
    for (i in seq_along(s11)) {
      x_pos <- max(df$rank.signature) * 1.1
      y_pos <- ypos[i]

      # text(x = x_pos, y = y_pos, labels = df$label2[s11[i]], pos = 4, cex = 0.7)
      text(x = x_pos, y = y_pos, labels = df$label2[s11[i]], pos = 4, cex = 0.7, family = "sans", font = 3)
      segments(x0 = df$rank.signature[s11[i]], y0 = df$Stouffer.signature[s11[i]], x1 = x_pos, y1 = y_pos, col = "lightgrey")
      segments(x0 = x_pos, y0 = y_pos, x1 = x_pos - 0.02 * max(df$rank.signature), y1 = y_pos, col = "lightgrey")
    }
  }


}

# Example usage:
# Assuming you have the variables 'signature', 'geneset', and 'gsea.out' already defined
#WaterfallGSEA2(signature = signature, geneset = geneset, gsea.out = gsea.out)
