
#' Install Python dependencies.
#' @export
install_py_dependencies <- function() {
    reticulate::py_install("pandas")
    reticulate::py_install("numpy")
    reticulate::py_install("plotly")
    reticulate::py_install("plotly_resampler")
}
