#' Generate a volcano plot of genes from WGCNA modules and a differential expression (limma) analysis
#'
#' Generate a volcano plot of genes from a differential expression (limma) analysis, with point
#' inclusion and color based on membership in modules from a WGCNA analysis. This plot can be output
#' to a plotting window, or to a pdf. Points can be selected for inclusion based on their connectivity
#' to the module eigengene (kME). Points can be labeled with gene names, and the points to be labeled
#' can be set based on an ellipse oriented to the x- and y-axes.
#' @param topGenes a data frame, typically the output of a call to \code{topTable}. Must contain genes, log2 fold-change, and adjusted p-values.
#' @param gene_modules a data frame containing the gene module memberships. Must contain, at minimum, columns for \code{gene}, \code{color_assigned}, and \code{kME_color_assigned}. If provided, the function outputs a pdf of the plot, named "{file_prefix}.pdf".
#' @param plotdims a numeric vector, the size (in inches) of the plotting object. Either the size of the pdf, or the size of the plotting window.
#' @param fc_cut numeric, the (absolute value) log2 fold-change threshold for determining significance of genes. This value is also plotted as vertical dotted lines.
#' @param p_cut numeric, the p-value threshold for determining significance of genes. This value is also plotted as a horizontal dotted line.
#' @param x_lim,y_lim numeric vectors, the lower and upper limits of the plotting space along the x- and y-axes. Passed to \code{ggplot2::xlim}.
#' @param gene_labs logical, whether to include gene labels for genes with extreme logFC and p-value. If \code{TRUE}, genes with values outside the labeling ellipse will be labeled.
#' @param x_cut,y_cut numeric, the radii of the labeling ellipse along the x- and y-axes. Genes with values outside the ellipse are labeled with gene names. Default to 0, which results in all genes being labeled.
#' @param point_order character string, specifying how to order the points. Currently accepted values are "random", which randomizes the order of the points, and "input", which send the points to ggplot as they are in the input data frame. Defaults to "random".
#' @import ggplot2
#' @import stringr
#' @details Genes with \code{kME_color_assigned > kME_cut} are plotted. Colors of points are determined by the values of \code{color_assigned}; if different colors are preferred, these values must be modified prior to calling \code{plot_WGCNA_module_volcano_2var}.
#' @export
#' @usage \code{
#' plot_WGCNA_module_volcano_2var(
#'      topGenes, gene_modules,
#'      file_prefix=NULL, plotdims=c(9,9),
#'      kME_cut=0.75, fc_cut=log2(1.5), p_cut=0.01,
#'      x_lim=NULL, y_lim=NULL,
#'      gene_labs=FALSE, x_cut=0, y_cut=0,
#'      point_order="random"
#'      )}
plot_WGCNA_module_volcano_2var <-
  function(topGenes, gene_modules,
           file_prefix=NULL, plotdims=c(9,9),
           kME_cut=0.75, fc_cut=log2(1.5), p_cut=0.01,
           x_lim=NULL, y_lim=NULL,
           gene_labs=FALSE, x_cut=0, y_cut=0,
           point_order="random") {
    topGenes <- topGenes[na.omit(match(gene_modules$gene, rownames(topGenes))),] # remove genes not in any module
    topGenes[,c("color_assigned", "kME_color_assigned")] <-
      gene_modules[match(rownames(topGenes), gene_modules$gene),
                   c("color_assigned", "kME_color_assigned")] # get colors and kME values for genes in WGCNA modules
    topGenes <- topGenes[topGenes$kME_color_assigned >= kME_cut,] # remove genes below kME connectivity to WGCNA module
    topGenes <- miscHelpers::order_points(topGenes, method=point_order)
    
    # generate volcano plot
    volcano <-
      ggplot(data = topGenes,
             aes(x=logFC, y=-log10(adj.P.Val),
                 colour=factor(color_assigned, levels=unique(color_assigned)))) +
      geom_point(alpha=0.6, size=3) +
      theme(legend.position = "none") +
      xlab("log2 fold change") + ylab("-log10 Adj P") +
      geom_vline(xintercept = fc_cut, linetype="dotted", size=1.0) +
    geom_vline(xintercept = -fc_cut, linetype="dotted", size=1.0) +
    geom_hline(yintercept = -log10(p_cut), linetype="dotted",size=1.0) + 
    scale_colour_manual(values=unique(topGenes$color_assigned))
  if (!is.null(x_lim)) {volcano <- volcano + xlim(x_lim)}
  if (!is.null(y_lim)) {volcano <- volcano + ylim(y_lim)}
  if (gene_labs) {
    volcano <- volcano +
      geom_text(
        data=topGenes[((topGenes$logFC^2)/(x_cut^2) +
                         (log10(topGenes$adj.P.Val)^2)/(y_cut^2)) > 1,],
        aes(label=rownames(topGenes[((topGenes$logFC^2)/(x_cut^2) +
                                       (log10(topGenes$adj.P.Val)^2)/(y_cut^2)) > 1,])),
        color="black", size=3, vjust=1, hjust=0.5)}
    
    # output volcano plot to file or plot window
  if (!is.null(file_prefix)) {
    pdf(file=paste(file_prefix, paste("WGCNA_colors", kME_cut, sep="_"), "pdf", sep="."),
        w=plotdims[1], h=plotdims[2])
    on.exit(dev.off(), add=TRUE) # close plotting device on exit
  } else quartz(plotdims[1],plotdims[2])
  
  print(volcano)
  
}
