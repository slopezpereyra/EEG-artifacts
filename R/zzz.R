iplot.eeg <- NULL
iplot.analysis <- NULL

.onLoad <- function(libname, pkgname) {
    question <- "Install Python dependencies? Downloads numpy,
                pandas, plotly and ploty_resampler if the packages are missing.
                This dependencies are necessary for interactive plotting."

    user_permission <- utils::askYesNo(question)

    if (isTRUE(user_permission)) {
        reticulate::py_install("pandas")
        reticulate::py_install("numpy")
        reticulate::py_install("plotly")
        reticulate::py_install("plotly_resampler")
    } else {
        message("To use interactive plotting, run
        `artifactor::install_py_dependencies()` first")
    }

    print("Setting IPLOTTER function")
    iplotter <- reticulate::import_from_path(module = "iplotter", path = system.file("python", package = packageName()))
    iplot.eeg <<- iplotter$plot_eeg
    iplot.analysis <<- iplotter$plot_analysis
}
