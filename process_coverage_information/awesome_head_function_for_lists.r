# awesome head function for e.g. data frames inside a list
# thanks to https://gist.github.com/pimentel/256fc8c9b5191da63819

head.list <- function(obj, n = 6L, ...)
{
    stopifnot(length(n) == 1L)
    origN <- n
    n <- if (n < 0L)
        max(length(obj) + n, 0L)
    else min(n, length(obj))
    lapply(obj[seq_len(n)], function(x)
           {
               tryCatch({
                   head(x, origN, ...)
               }, error = function(e) {
                   x
               })
           })
}
environment(head.list) <- asNamespace('utils')