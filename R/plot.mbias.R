#' Plot DAGs before and after conditioning on collider (M bias)
#'
#' Create two Directed Acyclic Graphs (DAGs), before and after conditioning on the
#' collider M, for selection bias caused by M bias, using 'ggdag'.
#'
#' @param x 'mbias' object to plot.
#' @param type DAG before or after conditioning on M.
#' @param dec Number of digits displayed.
#' @param ... Other unused arguments.
#'
#' @return A DAG for selection bias caused by M bias.
#'
#' @family visualization
#'
#' @seealso \code{\link{mbias}}
#'
#' @examples
#' plot(mbias(or = c(2, 5.4, 2.5, 1.5, 1),
#' var = c("HIV", "Circumcision", "Muslim", "Low CD4", "Participation")))
#'
#' @export
#' @importFrom ggplot2 ggplot aes fortify annotate theme theme_bw
#' @importFrom dagitty coordinates
#' @importFrom ggdag dagify coords2df coords2list m_bias confounder_triangle ggdag theme_dag
plot.mbias <- function(x,
                       type = c("before", "after"),
                       dec = 2,
                       ...) {
    obj <- x
    type <- match.arg(type)

    if (type == "before") {
        coords <- dagitty::coordinates(m_bias()) %>%
            coords2df()

        .dag <- dagify(x ~ a,
                       m ~ a + b,
                       y ~ x + b,
                       exposure = "x",
                       outcome = "y",
                       labels = c(x = obj[[5]][2],
                                  y = obj[[5]][1],
                                  a = obj[[5]][3],
                                  b = obj[[5]][4],
                                  m = obj[[5]][5]),
                       coords = coords2list(coords))

        ggdag(.dag, text = FALSE, use_labels = "label") +
            theme_dag() +
            annotate("text", x = .1, y = .5,
                     label = round(obj[[4]][1], dec)) +
            annotate("text", x = .4, y = .7,
                     label = round(obj[[4]][2], dec)) +
        annotate("text", x = 1.6, y = .7,
                 label = round(obj[[4]][3], dec)) +
        annotate("text", x = 1.9, y = .5,
                 label = round(obj[[4]][4], dec))
    }

    else {
        coords <- dagitty::coordinates(confounder_triangle()) %>%
            coords2df()

        .dag <- dagify(x ~ z,
                       y ~ x + z,
                       exposure = "x",
                       outcome = "y",
                       labels = c(x = obj[[5]][2],
                                  y = obj[[5]][1],
                                  z = obj[[5]][5]),
                       coords = coords2list(coords))

        ggdag(.dag, use_labels = "label") +
            theme_dag() +
            annotate("text", x = .4, y = .5,
                     label = round(obj[[2]][1], dec)) +
        annotate("text", x = 1.6, y = .5,
                     label = round(obj[[2]][2], dec)) +
        annotate("text", x = 1, y = 0.05,
                     label = round(obj[[2]][3], dec))
    }
}
