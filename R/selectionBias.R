selectionBias <-
function(case, exposed) {
    if(inherits(case, c("table", "matrix")))
        tab <- case
    else tab <- table(case, exposed)
    return(tab)
}
