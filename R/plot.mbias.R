#' Plot M bias
#'
#' @export
plot.mbias <- function(x,
                       title1 = "DAG before conditioning on P",
                       title2 = "DAG after conditioning on P",
                       title.size = 6,
                       size = 6,
                       dec = 2, 
                       layout = c("landscape", "portrait"),
                       ...) {
    layout <- match.arg(layout)
    
    res.df <- plyr::ldply(x[1:3], data.frame)
    labs <- x[[4]]

    library(grid)
    library(gridExtra)
    
    dag.before <- ggplot2::ggplot() +
        annotate("text", x = 3, y = Inf,
                 label = title1, vjust = 1.5, size = title.size) + 
        annotate("text", x = 5, y = 1, size = size, label = labs[1]) +
        annotate("text", x = 1, y = 1, size = size, label = labs[2]) +
        annotate("text", x = 1, y = 3, size = size, label = labs[3]) +
        annotate("text", x = 5, y = 3, size = size, label = labs[4]) +
        annotate("text", x = 3, y = 2, size = size, label = labs[5]) + 
        annotate("text", x = .5, y = 2, size = size,
                 label = round(res.df[5, 2], dec)) +
        annotate("text", x = 5.5, y = 2, size = size,
                 label = round(res.df[8, 2], dec)) + 
        annotate("text", x = 2.5, y = 2.5, size = size,
                 label = round(res.df[6, 2], dec)) +
        annotate("text", x = 3.5, y = 2.5, size = size,
                 label = round(res.df[7, 2], dec)) +
        annotate("text", x = 3, y = 1.1, size = size,
                 label = round(res.df[9, 2], dec)) +
        geom_segment(aes(x = 1, y = 2.9, xend = 1, yend = 1.1),
                     arrow = arrow(length = unit(0.25, "cm"))) +
        geom_segment(aes(x = 5, y = 2.9, xend = 5, yend = 1.1),
                     arrow = arrow(length = unit(0.25, "cm"))) +
        geom_segment(aes(x = 1, y = 2.9, xend = 2.75, yend = 2.1),
                     arrow = arrow(length = unit(0.25, "cm"))) +
        geom_segment(aes(x = 5, y = 2.9, xend = 3.25, yend = 2.1),
                     arrow = arrow(length = unit(0.25, "cm"))) +
        geom_segment(aes(x = 1.75, y = 1, xend = 4.25, yend = 1),
                     arrow = arrow(length = unit(0.25, "cm")),
                     linetype = "dashed") +
        xlim(0, 6) + ylim(.5, 3.5) + 
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              line = element_blank(),
              text = element_blank())

    dag.after <- ggplot2::ggplot() +
        annotate("text", x = 3, y = Inf,
                 label = title2, vjust = 1.5, size = title.size) + 
        annotate("text", x = 5, y = 1, size = size, label = labs[1]) +
        annotate("text", x = 1, y = 1, size = size, label = labs[2]) +
        annotate("text", x = 1, y = 3, size = size, label = labs[3]) +
        annotate("text", x = 5, y = 3, size = size, label = labs[4]) +
        annotate("text", x = 3, y = 2, size = size, label = labs[5]) + 
        annotate("text", x = 3, y = 1.1, size = size,
                     label = round(res.df[3, 2], dec)) +
        annotate("text", x = 2.5, y = 1.5, size = size,
                     label = round(res.df[1, 2], dec)) +
        annotate("text", x = 3.5, y = 1.5, size = size,
                     label = round(res.df[2, 2], dec)) +
        geom_segment(aes(x = 3, y = 1.9, xend = 1, yend = 1.1),
                     arrow = arrow(length = unit(0.25, "cm")),
                     linetype = "dashed") +
        geom_segment(aes(x = 3, y = 1.9, xend = 5, yend = 1.1),
                     arrow = arrow(length = unit(0.25, "cm")),
                     linetype = "dashed") +
        geom_segment(aes(x = 1.75, y = 1, xend = 4.25, yend = 1),
                     arrow = arrow(length = unit(0.25, "cm")),
                     linetype = "dashed") +
        xlim(0, 6) + ylim(.5, 3.5) + 
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              line = element_blank(),
              text = element_blank())

    if(layout == "landscape")
        grid.arrange(dag.before, dag.after, nrow = 1)
    else grid.arrange(dag.before, dag.after, nrow = 2)
}
