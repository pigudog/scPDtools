#!/usr/bin/env Rscript

exchangeable_loom_version <- "3.0.0"

isExchangeableLoom <- function(h5f) {
  attrs <- rhdf5::h5readAttributes(h5f, "/")
  version <- attrs[["LOOM_SPEC_VERSION"]]
  return(!is.null(version) && version == exchangeable_loom_version)
}

readDimNames <- function(h5f) {
  cell_attr <- rhdf5::h5readAttributes(h5f, "/col_attrs")[["CellID"]]
  gene_attr <- rhdf5::h5readAttributes(h5f, "/row_attrs")[["Gene"]]
  source <- rhdf5::h5readAttributes(h5f, "/")[["created_from"]]
  if (!is.null(source) && source == "anndata") {
    if (is.null(cell_attr)) {
      cell_attr <- "obs_names"
    }
    if (is.null(gene_attr)) {
      gene_attr <- "var_names"
    }
  }
  return(list(col = cell_attr, row = gene_attr))
}

readManifest <- function(h5f) {
  tryCatch(
    {
      manifest <- data.frame(t(rhdf5::h5read(h5f, "/attr/manifest")), stringsAsFactors = FALSE)
      colnames(manifest) <- c("loom_path", "dtype", "anndata_path", "sce_path")
      return(manifest)
    },
    error = function(e) {
      manifest <- NULL
      return(manifest)
    }
  )
}

nestedEnvAsList <- function(x) {
  out <- as.list(x)
  lapply(out, function(e) if (is.environment(e)) nestedEnvAsList(e) else e)
}

nestedEnv <- function(root_env, paths, v) {
  n <- length(paths)
  var <- paths[1]
  if (n == 1) {
    root_env[[var]] <- v
  } else {
    if (is.null(root_env[[var]])) {
      root_env[[var]] <- new.env(parent = emptyenv())
    }
    nestedEnv(root_env[[var]], paths[2:n], v)
  }
  invisible()
}

flattenNestedListToEnv <- function(x, e, prefix = NULL) {
  entry_names <- names(x)
  if (is.null(entry_names)) entry_names <- seq(length(x))
  sapply(entry_names, function(name) {
    if (is.null(prefix)) {
      full_name <- name
    } else {
      full_name <- paste(prefix, name, sep = "__")
    }
    if (is.list(x[[name]])) {
      flattenNestedListToEnv(x[[name]], e, full_name)
    } else {
      e[[full_name]] <- x[[name]]
    }
  })
  invisible()
}

makeManifest <- function(entries, dtype = "array", loom_prefix = "/attr/", anndata_prefix = "/uns/",
                         sce_prefix = "@metadata$") {
  n <- length(entries)
  if (is.list(entries)) {
    entry_names <- names(entries)
    is_scalar <- sapply(entries, function(x) is.vector(x) && length(x) == 1)
    dtypes <- ifelse(is_scalar, "scalar", "array")
  } else {
    entry_names <- entries
    dtypes <- rep(dtype, n)
  }
  loom_paths <- paste0(loom_prefix, entry_names)
  if (endsWith(loom_prefix, "[")) {
    loom_paths <- paste0(loom_paths, "]")
  }
  if (is.null(anndata_prefix)) {
    anndata_paths <- rep("", n)
  } else {
    if (anndata_prefix == "/obsm/X_") {
      ad_names <- tolower(entry_names)
    } else {
      ad_names <- gsub("__", "/", entry_names)
    }
    anndata_paths <- paste0(anndata_prefix, ad_names)
  }
  if (startsWith(sce_prefix, "@metadata$")) {
    entry_names <- gsub("__", "$", entry_names)
  }
  sce_paths <- paste0(sce_prefix, entry_names)
  return(data.frame(
    loom_path = loom_paths, dtype = dtypes, anndata_path = anndata_paths,
    sce_path = sce_paths, stringsAsFactors = FALSE
  ))
}

readExchangeableLoom <- function(filename, backed = TRUE, main_layer_name = NULL) {
  stopifnot(file.exists(filename), rhdf5::H5Fis_hdf5(filename))
  h5f <- rhdf5::H5Fopen(filename, flag = "H5F_ACC_RDONLY")
  if (!isExchangeableLoom(h5f)) {
    rhdf5::H5Fclose(h5f)
    return(LoomExperiment::import(filename, type = "SingleCellLoomExperiment"))
  }
  dim_names <- readDimNames(h5f)
  manifest <- readManifest(h5f)
  rhdf5::h5closeAll()

  suppressWarnings(scle <- LoomExperiment::import(
    filename,
    type = "SingleCellLoomExperiment",
    rownames_attr = dim_names$row,
    colnames_attr = dim_names$col
  ))
  if (!backed) {
    for (i in seq_along(SummarizedExperiment::assays(scle))) {
      SummarizedExperiment::assays(scle)[[i]] <- as(SummarizedExperiment::assays(scle)[[i]], "dgCMatrix")
    }
  }

  h5f <- rhdf5::H5Fopen(filename, flag = "H5F_ACC_RDONLY")

  # Add appropriate assay name
  mx_attrs <- rhdf5::h5readAttributes(h5f, "/matrix")
  if (!is.null(main_layer_name)) {
    names(SummarizedExperiment::assays(scle))[1] <- main_layer_name
  } else if ("assay" %in% names(mx_attrs)) {
    names(SummarizedExperiment::assays(scle))[1] <- mx_attrs["assay"]
  } else {
    names(SummarizedExperiment::assays(scle))[1] <- "counts"
  }

  if (!is.null(manifest)) {

    # Graphs are already handled by import(), just record entries
    is_graph <- (startsWith(manifest$loom_path, "/col_graphs/") |
      startsWith(manifest$loom_path, "/row_graphs/"))

    # Handle reducedDims
    is_rd <- startsWith(manifest$loom_path, "/attr/reducedDims__")
    rd_paths <- manifest$loom_path[is_rd]
    names(rd_paths) <- sub("^@reducedDims@listData[$]", "", manifest$sce_path[is_rd])
    SingleCellExperiment::reducedDims(scle) <- S4Vectors::SimpleList(lapply(rd_paths, function(path) {
      mat <- t(rhdf5::h5read(h5f, path))
      rownames(mat) <- colnames(scle)
      mat
    }))

    # Handle global attributes
    is_attr <- startsWith(manifest$loom_path, "/.attrs[")
    src_paths <- manifest$loom_path[is_attr]
    tgt_paths <- manifest$sce_path[is_attr]
    attr_names <- substr(src_paths, 9, nchar(src_paths) - 1)
    global_attrs <- rhdf5::h5readAttributes(h5f, "/")
    mtdt <- new.env(parent = emptyenv(), size = length(tgt_paths))
    for (i in seq_along(attr_names)) {
      v <- global_attrs[[attr_names[i]]]
      paths <- unlist(strsplit(sub("^@metadata[$]", "", tgt_paths[i]), "$", fixed = TRUE))
      nestedEnv(mtdt, paths, v)
      scle@metadata[[attr_names[i]]] <- NULL
    }

    # Handle global datasets
    is_ds <- (!(is_graph | is_rd | is_attr) & manifest$sce_path != "")
    src_paths <- manifest$loom_path[is_ds]
    tgt_paths <- manifest$sce_path[is_ds]
    for (i in seq_along(tgt_paths)) {
      v <- rhdf5::h5read(h5f, src_paths[i])
      paths <- unlist(strsplit(sub("^@metadata[$]", "", tgt_paths[i]), "$", fixed = TRUE))
      nestedEnv(mtdt, paths, v)
    }
    mtdt <- nestedEnvAsList(mtdt)
    for (name in names(mtdt)) {
      scle@metadata[[name]] <- mtdt[[name]]
    }
  }
  rhdf5::h5closeAll()

  return(as(scle, "SingleCellExperiment"))
}

writeExchangeableLoom <- function(sce, filename, main_layer = NULL, return_manifest = FALSE) {
  scle <- LoomExperiment::SingleCellLoomExperiment(sce)

  # Clean rowData and colData
  row_fct_idx <- sapply(SummarizedExperiment::rowData(scle), is.factor)
  SummarizedExperiment::rowData(scle)[row_fct_idx] <- lapply(
    SummarizedExperiment::rowData(scle)[row_fct_idx],
    function(x) type.convert(as.character(x), as.is = TRUE)
  )
  col_fct_idx <- sapply(SummarizedExperiment::colData(scle), is.factor)
  SummarizedExperiment::colData(scle)[col_fct_idx] <- lapply(
    SummarizedExperiment::colData(scle)[col_fct_idx],
    function(x) type.convert(as.character(x), as.is = TRUE)
  )

  # Handle reducedDims. Move embeddings out of reducedDims so they don't get
  # written to unwanted location by export().
  rdims <- SingleCellExperiment::reducedDims(scle)
  SingleCellExperiment::reducedDims(scle) <- S4Vectors::SimpleList()
  if (!S4Vectors::isEmpty(rdims)) {
    rdim_manifest <- makeManifest(
      names(rdims),
      dtype = "array",
      loom_prefix = "/attr/reducedDims__",
      anndata_prefix = "/obsm/X_",
      sce_prefix = "@reducedDims@listData$"
    )
  } else {
    rdim_manifest <- NULL
  }

  # Handle graphs. They get written by export() but we still need to record the paths.
  if (!S4Vectors::isEmpty(LoomExperiment::colGraphs(scle))) {
    colgraph_manifest <- makeManifest(
      names(LoomExperiment::colGraphs(scle)),
      dtype = "graph",
      loom_prefix = "/col_graphs/",
      anndata_prefix = "/uns/",
      sce_prefix = "@colGraphs$"
    )
  } else {
    colgraph_manifest <- NULL
  }
  if (!S4Vectors::isEmpty(LoomExperiment::rowGraphs(scle))) {
    rowgraph_manifest <- makeManifest(
      names(LoomExperiment::rowGraphs(scle)),
      dtype = "graph",
      loom_prefix = "/row_graphs/",
      anndata_prefix = NULL,
      sce_prefix = "@rowGraphs$"
    )
  } else {
    rowgraph_manifest <- NULL
  }

  # Handle metadata. Flatten nested lists to make export() happy. Scalars go
  # to /.attrs, others go to /attr
  if (length(S4Vectors::metadata(scle)) > 0) {
    mtdt <- new.env(parent = emptyenv())
    flattenNestedListToEnv(scle@metadata, mtdt)
    mtdt <- lapply(mtdt, function(x) {
      if (!is.numeric(x) && !is.character(x) && !is.logical(x)) {
        x <- type.convert(as.character(x), as.is = TRUE)
      }
      if ("class" %in% names(attributes(x))) {
        attributes(x)$class <- NULL
      }
      x
    })
    is_attr <- rep(FALSE, length(mtdt))
    names(is_attr) <- names(mtdt)
    for (i in seq_along(mtdt)) {
      x <- mtdt[[i]]
      if ((is.vector(x) || is.array(x)) && (length(x) == 1)) {
        is_attr[i] <- TRUE
        mtdt[[i]] <- x[1]
      }
    }
    is_attr <- is_attr & !grepl("__", names(mtdt))

    # Let export handle attributes
    S4Vectors::metadata(scle) <- mtdt[is_attr]

    excluded_from_manifest <- c(
      "LOOM_SPEC_VERSION", "CreationDate", "last_modified",
      "CreatedWith", "LoomExperiment-class", "created_from",
      "last_modified_by"
    )
    attr_names <- names(S4Vectors::metadata(scle))
    attr_names <- attr_names[!attr_names %in% excluded_from_manifest]
    if (length(attr_names) > 0) {
      attr_manifest <- makeManifest(
        attr_names,
        dtype = "scalar",
        loom_prefix = "/.attrs[",
        anndata_prefix = "/uns/",
        sce_prefix = "@metadata$"
      )
    } else {
      attr_manifest <- NULL
    }

    datasets <- mtdt[!is_attr]
    if (length(datasets) > 0) {
      dts_manifest <- makeManifest(
        datasets,
        dtype = NULL,
        loom_prefix = "/attr/",
        anndata_prefix = "/uns/",
        sce_prefix = "@metadata$"
      )
    } else {
      dts_manifest <- NULL
    }
  } else {
    attr_manifest <- NULL
    dts_manifest <- NULL
    datasets <- list()
  }

  manifest <- rbind(attr_manifest, dts_manifest, rdim_manifest, colgraph_manifest, rowgraph_manifest)

  # LoomExperiment currently handles sparse matrix incorrectly, so convert to dense for now
  for (assay_name in SummarizedExperiment::assayNames(scle)) {
    if (class(SummarizedExperiment::assays(scle)[[assay_name]]) == "dgCMatrix") {
      SummarizedExperiment::assays(scle)[[assay_name]] <- as.matrix(SummarizedExperiment::assays(scle)[[assay_name]])
    }
  }
  # Write to loom by LoomExperiment::export
  if (file.exists(filename)) file.remove(filename)
  suppressWarnings(LoomExperiment::export(
    scle,
    filename,
    matrix = ifelse(!is.null(main_layer) && main_layer %in% SummarizedExperiment::assayNames(scle),
      main_layer, SummarizedExperiment::assayNames(scle)[1]
    ),
    colnames_attr = "obs_names",
    rownames_attr = "var_names"
  ))
  rhdf5::h5closeAll()

  # Write extra bits
  h5f <- rhdf5::H5Fopen(filename)

  # Write extra global attributes
  rhdf5::h5writeAttribute(exchangeable_loom_version, h5f, "LOOM_SPEC_VERSION")
  rhdf5::h5writeAttribute("sce", h5f, "created_from")

  # Write column names of 'CellID' and 'Gene' as attributes of '/col_attrs' and '/row_attrs'
  h5g_ca <- rhdf5::H5Gopen(h5f, "/col_attrs")
  rhdf5::h5writeAttribute("obs_names", h5g_ca, "CellID")
  h5g_ra <- rhdf5::H5Gopen(h5f, "/row_attrs")
  rhdf5::h5writeAttribute("var_names", h5g_ra, "Gene")

  # Write primary asssay name as attribute of '/matrix'
  h5d_mx <- rhdf5::H5Dopen(h5f, "/matrix")
  rhdf5::h5writeAttribute("assay", h5d_mx, names(SummarizedExperiment::assays(scle))[1])

  # Write manifest
  rhdf5::h5createGroup(h5f, "/attr")
  if (!is.null(manifest)) {
    manifest <- manifest[order(manifest$dtype, manifest$loom_path), ]
    rhdf5::h5write(t(manifest), h5f, "/attr/manifest")
  }

  # Write reducedDims
  for (i in seq_along(rdims)) {
    rdim <- rdims[[i]]
    loom_path <- rdim_manifest$loom_path[i]
    rhdf5::h5write(t(rdim), h5f, loom_path)
  }

  # Write extra global datasets
  for (i in seq_along(datasets)) {
    dts <- datasets[[i]]
    loom_path <- dts_manifest$loom_path[i]
    rhdf5::h5write(dts, h5f, loom_path)
  }

  # Remove '/col_attrs/reducedDims' to make anndata happy
  rhdf5::h5delete(h5f, "/col_attrs/reducedDims")
  rhdf5::h5closeAll()

  if (return_manifest) {
    return(manifest)
  } else {
    invisible()
  }
}
