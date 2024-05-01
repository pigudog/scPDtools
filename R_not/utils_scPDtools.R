

#' Download File from the Internet
#'
#' @inheritParams utils::download.file
#' @param methods Methods to be used for downloading files. The default is to try different download methods in turn until the download is successfully completed.
#' @param max_tries Number of tries for each download method.
#' @param ... Other arguments passed to \code{\link[utils]{download.file}}
#'
#' @importFrom utils download.file
#' @export
download <- function(url, destfile, methods = c("auto", "wget", "libcurl", "curl", "wininet", "internal"), quiet = FALSE, ..., max_tries = 2) {
  if (missing(url) || missing(destfile)) {
    stop("'url' and 'destfile' must be both provided.")
  }
  ntry <- 0
  status <- NULL
  while (is.null(status)) {
    for (method in methods) {
      status <- tryCatch(expr = {
        suppressWarnings(download.file(url, destfile = destfile, method = method, quiet = quiet, ...))
        status <- 1
      }, error = function(error) {
        message(error)
        message("Cannot download from the url: ", url)
        message("Failed to download using \"", method, "\". Retry...\n")
        Sys.sleep(1)
        return(NULL)
      })
      if (!is.null(status)) {
        break
      }
    }
    ntry <- ntry + 1
    if (is.null(status) && ntry >= max_tries) {
      stop("Download failed.")
    }
  }
  return(invisible(NULL))
}

#' Check and install R packages
#'
#' @param packages Package to be installed. Package source can be CRAN, Bioconductor or Github, e.g. scmap, davidsjoberg/ggsankey.
#' @param package_names The name of the package that corresponds to the \code{packages} parameter, used to check if the package is already installed.
#' By default, the package name is extracted according to the \code{packages} parameter.
#' @param install_methods Functions used to install R packages.
#' @param lib  The location of the library directories where to install the packages.
#' @param force Whether to force the installation of packages. Default is \code{FALSE}.
#'
#' @importFrom rlang %||%
#' @importFrom utils packageVersion
#' @export
check_R <- function(packages, package_names = NULL, install_methods = c("BiocManager::install", "install.packages", "devtools::install_github"), lib = .libPaths()[1], force = FALSE) {
  if (length(package_names) != 0 && length(package_names) != length(packages)) {
    stop("package_names must be NULL or a vector of the same length with packages")
  }
  status_list <- list()
  for (n in seq_along(packages)) {
    pkg <- packages[n]
    pkg_info <- pkg
    if (!grepl("/", pkg_info)) {
      pkg_info <- paste0("/", pkg_info)
    }
    if (!grepl("@", pkg_info)) {
      pkg_info <- paste0(pkg_info, "@")
    }
    git <- grep("/", sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\1", x = pkg_info), value = TRUE)
    git <- gsub("/", "", git)
    pkg_name <- package_names[n] %||% sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\2", x = pkg_info)
    version <- grep("@", sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\3", x = pkg_info), value = TRUE)
    version <- gsub("@", "", version)
    if (version != "") {
      force_update <- isTRUE(packageVersion(pkg_name) < package_version(version)) || isTRUE(force)
    } else {
      force_update <- isTRUE(force)
    }
    if (!suppressPackageStartupMessages(requireNamespace(pkg_name, quietly = TRUE)) || isTRUE(force_update)) {
      message("Install package: \"", pkg_name, "\" ...")
      status_list[[pkg]] <- FALSE
      i <- 1
      while (isFALSE(status_list[[pkg]])) {
        tryCatch(expr = {
          if (grepl("BiocManager", install_methods[i])) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager", lib = lib)
            }
            eval(str2lang(paste0(install_methods[i], "(\"", pkg, "\", lib=\"", lib, "\", update = FALSE, upgrade = \"never\", ask = FALSE, force = TRUE)")))
          } else if (grepl("devtools", install_methods[i])) {
            if (!requireNamespace("devtools", quietly = TRUE)) {
              install.packages("devtools", lib = lib)
            }
            if (!requireNamespace("withr", quietly = TRUE)) {
              install.packages("withr", lib = lib)
            }
            eval(str2lang(paste0("withr::with_libpaths(new = \"", lib, "\", ", install_methods[i], "(\"", pkg, "\", upgrade = \"never\", force = TRUE))")))
          } else {
            eval(str2lang(paste0(install_methods[i], "(\"", pkg, "\", lib=\"", lib, "\", force = TRUE)")))
          }
        }, error = function(e) {
          status_list[[pkg]] <- FALSE
        })
        if (version == "") {
          status_list[[pkg]] <- requireNamespace(pkg_name, quietly = TRUE)
        } else {
          if (requireNamespace(pkg_name, quietly = TRUE)) {
            status_list[[pkg]] <- packageVersion(pkg_name) >= package_version(version)
          } else {
            status_list[[pkg]] <- FALSE
          }
        }
        i <- i + 1
        if (i > length(install_methods)) {
          break
        }
      }
    } else {
      status_list[[pkg]] <- TRUE
    }
  }
  out <- sapply(status_list, isTRUE)
  out <- out[!out]
  if (length(out) > 0) {
    stop("Failed to install the package(s): ", paste0(names(out), collapse = ","), ". Please install manually.")
  }
}

kegg_get <- function(url) {
  temp <- tempfile()
  on.exit(unlink(temp))
  download(url = url, destfile = temp)
  content <- as.data.frame(do.call(rbind, strsplit(readLines(temp), split = "\t")))
  return(content)
}

rescale <- function(x, from = range(x, na.rm = TRUE, finite = TRUE), to = c(0, 1)) {
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  } else {
    return((x - from[1]) / diff(from) * diff(to) + to[1])
  }
}

zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
  if (length(x) == 1) {
    return(TRUE)
  }
  if (length(x) != 2) {
    stop("x must be length 1 or 2")
  }
  if (any(is.na(x))) {
    return(NA)
  }
  if (x[1] == x[2]) {
    return(TRUE)
  }
  if (all(is.infinite(x))) {
    return(FALSE)
  }
  m <- min(abs(x))
  if (m == 0) {
    return(FALSE)
  }
  abs((x[1] - x[2]) / m) < tol
}

#' @importFrom grDevices col2rgb rgb
col2hex <- function(cname) {
  colMat <- col2rgb(cname)
  rgb(red = colMat[1, ] / 255, green = colMat[2, ] / 255, blue = colMat[3, ] / 255)
}

#' Invoke a function with a list of arguments
#' @param .fn A function, or function name as a string.
#' @param .args A list of arguments.
#' @param Other arguments passed to the function.
#' @param .env Environment in which to evaluate the call. This will be most useful if .fn is a string, or the function has side-effects.
#' @importFrom rlang caller_env is_null is_scalar_character is_character is_function set_names env env_get env_bind syms call2
#' @export
invoke <- function(.fn, .args = list(), ..., .env = caller_env()) {
  args <- c(.args, list(...))
  .bury <- c(".fn", "")
  if (is_null(.bury) || !length(args)) {
    if (is_scalar_character(.fn)) {
      .fn <- env_get(.env, .fn, inherit = TRUE)
    }
    call <- call2(.fn, !!!args)
    return(.External2(rlang:::ffi_eval, call, .env))
  }
  if (!is_character(.bury, 2L)) {
    abort("`.bury` must be a character vector of length 2")
  }
  arg_prefix <- .bury[[2]]
  fn_nm <- .bury[[1]]
  buried_nms <- paste0(arg_prefix, seq_along(args))
  buried_args <- set_names(args, buried_nms)
  .env <- env(.env, !!!buried_args)
  args <- set_names(buried_nms, names(args))
  args <- syms(args)
  if (is_function(.fn)) {
    env_bind(.env, `:=`(!!fn_nm, .fn))
    .fn <- fn_nm
  }
  call <- call2(.fn, !!!args)
  .External2(rlang:::ffi_eval, call, .env)
}

#' Implement similar functions to the \code{unnest} function in the tidyr package
#' @param data A data frame.
#' @param cols Columns to unnest.
#' @param keep_empty By default, you get one row of output for each element of the list your unchopping/unnesting. This means that if there's a size-0 element (like \code{NULL} or an empty data frame), that entire row will be dropped from the output. If you want to preserve all rows, use \code{keep_empty = TRUE} to replace size-0 elements with a single row of missing values.
#' @export
unnest <- function(data, cols, keep_empty = FALSE) {
  if (nrow(data) == 0 || length(cols) == 0) {
    return(data)
  }
  for (col in cols) {
    col_expand <- unlist(data[[col]])
    expand_times <- sapply(data[[col]], length)
    if (isTRUE(keep_empty)) {
      data[[col]][expand_times == 0] <- NA
      col_expand <- unlist(data[[col]])
      expand_times[expand_times == 0] <- 1
    }
    data <- data[rep(seq_len(nrow(data)), times = expand_times), ]
    data[, col] <- col_expand
  }
  rownames(data) <- NULL
  return(data)
}

#' Attempts to turn a dgCMatrix into a dense matrix
#' @param matrix A dgCMatrix
#' @export
as_matrix <- function(matrix) {
  if (!inherits(matrix, "dgCMatrix")) {
    stop("matrix is not a dgCMatrix.")
  }
  row_pos <- matrix@i
  col_pos <- findInterval(seq(matrix@x) - 1, matrix@p[-1])

  out <- asMatrix(
    rp = row_pos, cp = col_pos, z = matrix@x,
    nrows = matrix@Dim[1], ncols = matrix@Dim[2]
  )

  row.names(out) <- matrix@Dimnames[[1]]
  colnames(out) <- matrix@Dimnames[[2]]
  return(out)
}

#' Capitalizes the characters
#' Making the first letter uppercase
#'
#' @examples
#' x <- c("dna methylation", "rRNA processing", "post-Transcriptional gene silencing")
#' capitalize(x)
#' @param x A vector of character strings to be capitalized.
#' @param force_tolower Whether to force the remaining letters to be lowercase.
#' @export
capitalize <- function(x, force_tolower = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  if (!inherits(x, "character")) {
    stop("x must be the type of character.")
  }
  if (isTRUE(force_tolower)) {
    x <- paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))), sep = "")
  } else {
    first_word <- sapply(strsplit(x, "\\s|-"), function(s) s[1])
    index <- which(first_word == tolower(first_word))
    x[index] <- paste(toupper(substr(x[index], 1, 1)), substr(x[index], 2, nchar(x[index])), sep = "")
  }
  return(x)
}

str_wrap <- function(x, width = 80) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  x_wrap <- unlist(lapply(x, function(i) paste0(strwrap(i, width = width), collapse = "\n")))
  return(x_wrap)
}

#' Split a vector into the chunks
#'
#' @param x A vector.
#' @param nchunks Number of chunks.
#' @examples
#' x <- 1:10
#' names(x) <- letters[1:10]
#' tochunks(x, nchunks = 3)
#' @export
tochunks <- function(x, nchunks) {
  split(x, cut(seq_along(x), nchunks, labels = FALSE))
}

#' Generate a iterator along chunks of a vector
#' @param x A vector.
#' @param nchunks Number of chunks.
#' @examples
#' \dontrun{
#' library(BiocParallel)
#' x <- 1:100
#' BPPARAM <- bpparam()
#' bpprogressbar(BPPARAM) <- TRUE
#' bpworkers(BPPARAM) <- 10
#' slow_fun <- function(x) {
#'   out <- NULL
#'   for (i in seq_along(x)) {
#'     Sys.sleep(0.5)
#'     out[[i]] <- x[[i]] + 3
#'   }
#'   return(out)
#' }
#' system.time({
#'   res0 <- lapply(x, FUN = slow_fun)
#' })
#' unlist(res0, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res1 <- bplapply(x, FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res1, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res2 <- bplapply(tochunks(x, nchunks = bpworkers(BPPARAM)), FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res2, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res3 <- bpiterate(ITER = iterchunks(x, nchunks = bpworkers(BPPARAM)), FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res3, recursive = FALSE, use.names = FALSE)[71:73]
#' }
#' @export
iterchunks <- function(x, nchunks) {
  chunks <- tochunks(x, nchunks)
  i <- 0L
  function() {
    if (i >= length(chunks)) {
      return(NULL)
    }
    i <<- i + 1L
    x[chunks[[i]]]
  }
}
