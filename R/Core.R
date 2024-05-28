# Main function and helper functions for scgft


#### Main function ####

#' @title Generates synthetic scRNA-seq data in-silico
#' @description   Leverages the principles of Fourier Transform, to generate
#'   synthetic cells, while preserving the biological plausibility of the
#'   original population.
#'
#' @param    object            A seurat object.
#' @param    nsynth            Number of synthetic cells to be generated.
#' @param    ncpmnts           The number of components to modify during the
#'                             synthesis process.
#' @param    groups            A character in the original \code{object} metadata used
#'                             for group-based synthesis (such as cell types or clusters).
#'                             All groups are uniformly scaled for desired expansion.
#' @param    cells             Specifies the barcode(s) of the cell(s) to be
#'                             used for cell-based synthesis. If a \code{list} of
#'                             barcodes is provided, \code{nsynth} cells will be
#'                             synthesized for each barcode in the list. If a
#'                             \code{vector} of barcodes is provided, \code{nsynth}
#'                             cells will be synthesized for the specified group of barcodes.
#'
#' @return Returns a combined \code{object} of original and synthetic cells. The
#'   synthesized cells are integrated into the \code{data} layer, while their
#'   adjusted counts are integrated into the \code{counts} layer.
#'   The synthesized cells are labeled by the cell barcodes they originated
#'   from, suffixed with "_synth". A new \code{synthesized} column is added to
#'   metadata indicating the synthesis status of each cell, labeled as either
#'   "yes" or "no".
#'
#' @details The input \code{object} should be log normalized and contain top
#'   variable features.
#'
#' @export
RunScGFT <- function(object, nsynth, ncpmnts=1,
                     groups=NULL, cells=NULL) {
  # =======================================
  # Check if Seurat package is installed
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running scGFT requires 'Seurat'.")
  }
  if (!requireNamespace('SeuratObject', quietly = TRUE)) {
    stop("Running scGFT requires 'SeuratObject'.")
  }
  # Check if Seurat object contains RNA assay
  if (!"RNA" %in% Seurat::Assays(object)) {
    stop("Object does not contain RNA assay.")
  }
  # Check if Seurat object contains RNA assay with layers for counts and data
  if (!all(c("counts", "data") %in% SeuratObject::Layers(object, assay="RNA"))) {
    stop("Layers for counts or data does not exist. Please run Seurat::NormalizeData() first.")
  }
  # Check if NN graphs is available
  if (!"RNA_nn" %in% Graphs(object)) {
    stop("Object should contain 'RNA_nn' graph. Please use `Seurat::FindNeighbors()'.")
  }
  # Check if Seurat object contains variable features
  if (is.null(Seurat::VariableFeatures(object))) {
    stop("Seurat object does not contain variable features. Please run Seurat::FindVariableFeatures() first.")
  }
  # check if 'groups' exist in the metadata
  if (is.null(cells) & is.null(groups)) {
    stop(paste("One of 'groups' or 'cells'", "should be assigned."))
  }
  # check if 'groups' exist in the metadata
  if (is.null(cells) & !is.null(groups)) {
    if (!groups %in% names(object@meta.data)) {
      stop(paste0("'", groups, "' does not exist in the object metadata."))
    }
  }
  # check ncpmnts and number of variable features
  nvf <- length(Seurat::VariableFeatures(object))
  if (ncpmnts > nvf) {
    stop(paste0("'ncpmnts' (", ncpmnts, ") cannot be larger than number of variable features (", nvf, ")."))
  }
  # =======================================
  suppressWarnings({
    orig_cnt <- as.matrix(Seurat::GetAssayData(object, assay="RNA", layer="counts"))
    })
  suppressWarnings({
    orig_dta <- as.matrix(Seurat::GetAssayData(object, assay="RNA", layer="data"))
    })
  suppressWarnings({
    adj_mtx <- as.matrix(attr(object, which="graphs")[["RNA_nn"]])
  })
  # =======================================
  genes <- rownames(object)
  varftrs <- Seurat::VariableFeatures(object)
  invarftrs <- setdiff(genes, varftrs)
  # =======================================
  start_time <- Sys.time()
  # =======================================
  if (!is.null(cells)) {
    if (is.atomic(cells)) {
      var_mtx <- orig_dta[varftrs, cells, drop=FALSE]
      invar_mtx <- orig_dta[invarftrs, cells, drop=FALSE]
      groups <- rep(1, length(cells))
      syn_mtx <- PerformDIFT(var_mtx=var_mtx, groups=groups, nsynth=nsynth, ncpmnts=ncpmnts, adj_mtx=adj_mtx)
    } else if (is.list(cells)) {
      syn_mtx <- lapply(cells, function(x) {
        var_mtx <- orig_dta[varftrs, x, drop=FALSE]
        groups <- rep(1, length(x))
        PerformDIFT(var_mtx=var_mtx, groups=groups, nsynth=nsynth, ncpmnts=ncpmnts, adj_mtx=adj_mtx)
      })
      syn_mtx <- do.call(cbind, syn_mtx)
      invar_mtx <- orig_dta[invarftrs, unlist(cells), drop=FALSE]
    } else {
      stop()
    }
  } else if (!is.null(groups)){
    var_mtx <- orig_dta[varftrs, ]
    invar_mtx <- orig_dta[invarftrs, ]
    groups <- as.character(object@meta.data[[groups]])
    syn_mtx <- PerformDIFT(var_mtx=var_mtx, groups=groups, nsynth=nsynth, ncpmnts=ncpmnts, adj_mtx=adj_mtx)
  } else {
    stop()
  }
  # =======================================
  end_time <- Sys.time()
  dt <- as.numeric(difftime(end_time, start_time, units = "secs"))
  message(paste("Synthesis completed in:", round(dt/60, 2), "min"))
  # ===================================
  # get the synthesized cells names
  synt_cells_nm <- FetchSynthesizedCells(colnames(syn_mtx))
  # ===================================
  message(paste("Integrating data (1/2)"))
  # extend to integrated matrix
  syn_mtx_full <- ExtendSynthesizedMatrix(syn_mtx, genes, invar_mtx, invarftrs, synt_cells_nm)
  metadata_synt <- ExtendMetaData(object@meta.data, syn_mtx_full)
  # ===================================
  message(paste("Integrating data (2/2)"))
  sobj_synt <- SobjMerger(cnt_ls = list(as(orig_dta, "dgCMatrix"), as(syn_mtx_full, "dgCMatrix")),
                          mtd_ls = list(object@meta.data, metadata_synt))
  sobj_synt@assays$RNA$data <- sobj_synt@assays$RNA$counts
  sobj_synt@assays$RNA$counts <- GetCountMatrix(sobj_synt, orig_cnt, orig_dta, synt_cells_nm, syn_mtx_full)
  # ===================================
  Seurat::VariableFeatures(sobj_synt) <- varftrs
  # ===================================
  ids <- grepl("_synth", rownames(sobj_synt@meta.data))
  sobj_synt@meta.data$synthesized <- "no"
  sobj_synt@meta.data$synthesized[ids] <- "yes"
  # ===================================
  x <- format(dim(sobj_synt)[2], big.mark=",")
  y <- format(sum(sobj_synt$synthesized == "yes"), big.mark=",")
  message(paste("A Seurat object with", x, "cells, including", y, "synthesized."))
  # ===================================
  return(sobj_synt)
}


#### Summary functions ####


#' @title synthesized cells QC
#' @description   Calculates the likelihood that synthesized cells will have the
#'   same identity as their original counterparts.
#'
#' @param    object   A Seurat object that includes synthesized cells, i.e., an object post \link{RunScGFT}.
#' @param    groups   Same character in the original \code{object} metadata used for synthesis by \link{RunScGFT}.
#'
#' @return   synthesis accuracy and deviation from original cells in percentage.
#'
#' @details assesses the quality of the synthesis process by calculating the
#'   proportion of synthesized cells that correctly grouped with their
#'   corresponding original cells according to the given \code{groups} used for synthesis.
#'
#' @export
#'
statsScGFT <- function(object, groups) {
  # ===================================
  # Hack for visibility of dplyr variables
  . <- NULL
  barcode <- orig_cell <- synthesized <- NULL
  original_barcode <- original_id <- NULL
  synthesized_barcode <- synthesized_id <- NULL
  cluster_match <- matching_clusters <- total_cells <- accuracy_percentage <- NULL
  all_of <- NULL
  # ======================
  # check if 'groups' exist in the metadata
  if (!groups %in% names(object@meta.data)) {
    stop(paste0("'", groups, "' does not exist in the object metadata."))
  }
  if (!"synthesized" %in% names(object@meta.data)) {
    stop(paste("'synthesized' does not exist in the object metadata."))
  }
  # ===================================
  vars <- c("barcode", groups, "orig_cell", "synthesized")
  # ===================================
  stats_df <- object@meta.data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::mutate(orig_cell = sub("_synth.*$", "", barcode))
  # ===================================
  originals <- stats_df[, colnames(stats_df) %in% vars] %>%
    dplyr::filter(synthesized == "no") %>%
    dplyr::rename(original_barcode = barcode,
                  original_id = groups)
  # ===================================
  synths <- stats_df[, colnames(stats_df) %in% vars] %>%
    dplyr::filter(synthesized == "yes") %>%
    dplyr::rename(synthesized_barcode = barcode,
                  synthesized_id = groups)
  # ===================================
  results_stats <- merge(synths, originals, by = "orig_cell", all.x = TRUE) %>%
    dplyr::mutate(cluster_match = original_id == synthesized_id) %>%
    dplyr::summarise(matching_clusters = sum(cluster_match),
                     total_cells = n(),
                     accuracy_percentage = (matching_clusters / total_cells) * 100)
  # ===================================
  acr <- round(results_stats$accuracy_percentage, 2)
  message(paste("Synthesized cells:", format(results_stats$total_cells, big.mark=",")))
  message(paste("Matching groups:", format(results_stats$matching_clusters, big.mark=",")))
  message(paste("Accuracy (%):", acr))
  # ===================================
  # message("Calculating deviation from originals...")
  # devs <- DevScGFT(object)
  # message(paste("Deviation (%):", round(mean(devs)*100, 2), "+/-", round(sd(devs)*100, 2)))
  # ===================================
  return(list("accuracy" = acr))
}


#### Helper functions ####


# Initialize progress bar function
initialize_bar <- function(totbar, wdth=100) {
  progress_bar$new(
    format = "  [:bar] :percent in :elapsed",
    total = totbar, clear = FALSE, width = wdth
  )
}


# Merge seurat objects
SobjMerger <- function(cnt_ls, mtd_ls) {
  # ===================================
  mtd_ls <- lapply(mtd_ls, function(x){
    x[] <- lapply(x, function(y) if (is.factor(y)) as.character(y) else y)
    return(x)
  })
  mtd <- Reduce(rbind, mtd_ls)
  # ===================================
  # message(paste("*** create merged seurat object..."))
  sobj <- Seurat::CreateSeuratObject(counts=cnt_ls, meta.data=mtd)
  # ===================================
  # message(paste("*** join all layers..."))
  sobj <- SeuratObject::JoinLayers(sobj, assay="RNA")
  # ===================================
  return(sobj)
}


# Sample a modification factor from a Gaussian distribution
# with mean and sd for the component across the group
# ModificationGroup <- function(ft_mtx, cmpnts, cells_in_grp) {
#   ampltds <- abs(ft_mtx[cmpnts, cells_in_grp, drop = FALSE])
#   ampltds <- apply(ampltds, 1, function(x) x/sqrt(sum(x^2)))
#   mdfs <- apply(ampltds, 2, function(x) rnorm(1, mean=2, sd=sd(x))) # the previous line flips the dimensions
#   # stopifnot(sum(is.na(mdfs)) == 0)
#   return(mdfs)
# }


ModificationGroup <- function(cmpnts, distr_ls, grp_id) {
  sapply(distr_ls[[grp_id]]$sds[cmpnts], function(x){
    rnorm(1, mean=2, sd=x)
  })
}


# Sample a modification factor from a Normal Gaussian distribution
ModificationSingle <- function(cmpnts) {
  rnorm(length(cmpnts), mean=2, sd=1)
}


# get modes of variation present in groups of cells
modesVar <- function(sorted_grp_ids, cells_bc, groups, ft_mtx) {
  distr_ls <- lapply(sorted_grp_ids, function(grp_id){
    cells_in_grp <- cells_bc[groups == grp_id]
    ampltds <- abs(ft_mtx[, cells_in_grp, drop=FALSE])
    ampltds <- apply(ampltds, 1, function(x) x/sqrt(sum(x^2)))
    list(mns = apply(ampltds, 2, mean), # the previous line flips the dimensions,
         sds = apply(ampltds, 2, sd)) # the previous line flips the dimensions)
  })
  names(distr_ls) <- sorted_grp_ids
  return(distr_ls)
}


# Performs Discrete and INverse Fourier Transforms
PerformDIFT <- function(var_mtx, groups, nsynth, ncpmnts=1, adj_mtx) {
  # Extracting components from the input object
  genes_nm <- rownames(var_mtx)
  cells_bc <- colnames(var_mtx)
  # originalUMICount <- colSums(var_mtx) # Obtain original cell's total UMI count for normalization
  # =======================================
  # Step 1: Calculate the number of cells to synthesize per group
  groups_cnt <- table(groups)
  nper_grp <- round(nsynth*groups_cnt/sum(groups_cnt))
  # =======================================
  excess_cells <- sum(nper_grp) - nsynth
  # Adjust the group with the maximum number of cells if there's an excess
  if (excess_cells > 0) {
    max_group <- which.max(nper_grp)
    nper_grp[max_group] <- nper_grp[max_group] - excess_cells
  }
  # Similarly, you can adjust for a deficit in the total (if sum(nper_grp) < nsynth)
  # by adding missing cells to the group with the least cells
  missing_cells <- nsynth - sum(nper_grp)
  if (missing_cells > 0) {
    min_group <- which.min(nper_grp)
    nper_grp[min_group] <- nper_grp[min_group] + missing_cells
  }
  # =======================================
  # geneMatrix: Matrix with gene names as rows and cell barcodes as columns
  # groups: Vector containing groups IDs for each cell
  # Step 1: Apply Fourier Transform to each cell
  message(paste("Discrete fourier transform..."))
  ft_mtx <- apply(var_mtx, 2, function(cell) fft(cell))
  rownames(ft_mtx) <- paste0("cmpnt_",1:nrow(var_mtx))
  len_ft <- nrow(ft_mtx)
  # =======================================
  sorted_grp_ids <- as.character(names(sort(nper_grp, decreasing=T)))
  # =======================================
  # get modes of variation present in groups of cells
  if (dim(ft_mtx)[2] > 1) {
    distr_ls <- modesVar(sorted_grp_ids, cells_bc, groups, ft_mtx)
  }
  # =======================================
  message(paste("Inverse fourier transform..."))
  message(paste("synthesizing", format(sum(nper_grp), big.mark=","), "cells..."))
  # =======================================
  # Step 2: Generate synthetic cells for each group
  # set.seed(12345)
  synt_ls <- list()
  devs_c <- c()
  cell_cnt <- 0
  for (grp_id in sorted_grp_ids) {
    # get cells in the group
    cells_in_grp <- cells_bc[groups == grp_id]
    # decide modification function
    if (length(cells_in_grp) > 1) {
      ModificationFunc <- function(cmpnts) ModificationGroup(cmpnts, distr_ls, grp_id)
    } else {
      ModificationFunc <- function(cmpnts) ModificationSingle(cmpnts)
    }
    # Randomly select a cell from the group
    n_synt <- nper_grp[names(nper_grp) == grp_id]
    # cells_smpld <- sampleCells(n_synt, cells=cells_in_grp)
    cells_smpld <- sample(cells_in_grp, n_synt, replace=TRUE)
    synt_res <- lapply(cells_smpld, function(x){
      # get the selected cell FT
      modifiedFT <- ft_mtx[, x]
      # Randomly select a component to modify
      cmpnts <- sample(2:len_ft, ncpmnts, replace=FALSE) # do not modify Direct Current (DC) Component
      # Calculate the modification factor
      modFactors <- ModificationFunc(cmpnts)
      # Apply modification to the selected component
      modifiedFT[cmpnts] <- modifiedFT[cmpnts] * modFactors
      # Apply modification to the selected component conjugate
      boolid <- cmpnts != 1 & (len_ft %% 2 == 0 & cmpnts != len_ft/2 + 1) # here is a vector of boolens so used & instead of &&
      conjugateComps <- len_ft - cmpnts[boolid] + 2
      modifiedFT[conjugateComps] <- modifiedFT[conjugateComps] * modFactors[boolid]
      # Perform IFT to synthesize a new cell
      syn_cell <- Re(fft(modifiedFT, inverse=TRUE) / len_ft)
      # calculate deviations
      denom <- var_mtx[, x]
      denom[denom == 0] <- 1
      devs <- mean((pmax(0, syn_cell) - var_mtx[, x])/denom)
      # Find the indices of the neighbors. Combine the synthesized cell, its original
      # counterpart, and its neighbors. Calculate the average expression per gene.
      # rowMeans(cbind(syn_cell, var_mtx[, intersect(names(which(adj_mtx[x, ] > 0)), cells_in_grp), drop=FALSE]))
      syn_cell <- rowMeans(cbind(syn_cell, var_mtx[, names(which(adj_mtx[x, cells_in_grp] > 0)), drop=FALSE]))
      # return
      list("synt" = syn_cell,
           "devs" = devs)
    })
    # =======================================
    # extract cells
    synt <- lapply(synt_res, function(x) x$synt)
    names(synt) <- cells_smpld
    synt_ls <- c(synt_ls, synt)
    # =======================================
    # extract devs
    devs <- sapply(synt_res, function(x) x$devs)
    devs_c <- c(devs_c, devs)
    # =======================================
    cell_cnt <- cell_cnt + length(synt)
    message(paste(format(cell_cnt, big.mark=","), "cells synthesized..."))
  }
  stopifnot(length(synt_ls) == nsynth)
  stopifnot(unique(lengths(synt_ls)) == length(genes_nm))
  syn_mtx <- convert_list_to_matrix(synt_ls)
  rownames(syn_mtx) <- genes_nm
  syn_mtx[syn_mtx < 1e-6] <- 0 # Apply ReLU
  # =======================================
  message(paste("Deviation from originals (%):", round(mean(devs_c)*100, 2), "+/-", round(sd(devs_c)*100, 2)))
  # =======================================
  return(syn_mtx)
}


# Generated the raw count matrix by reversing log-transformed data
GetCountMatrix <- function(sobj_synt, orig_cnt, orig_dta, synt_cells_nm, syn_mtx_full) {
  # ===================================
  cnt_full <- matrix(NA, nrow=nrow(sobj_synt), ncol=ncol(sobj_synt))
  rownames(cnt_full) <- rownames(sobj_synt)
  colnames(cnt_full) <- colnames(sobj_synt)
  # ===================================
  indcs <- match(colnames(orig_cnt), colnames(cnt_full))
  cnt_full[, indcs] <- orig_cnt
  # ===================================
  # Initialize progress bar
  nl <- length(synt_cells_nm)
  pb <- initialize_bar(totbar=nl, wdth=66)
  # ===================================
  # Add count matrix to the full synth matrix, ensuring correct cell name alignment
  lents <- lengths(synt_cells_nm)
  for (x in names(synt_cells_nm)) {
    pb$tick() # Update progress bar
    cls <- synt_cells_nm[[x]]
    denom <- y <- replicate(lents[x], orig_dta[, x])
    denom[denom == 0] <- 1
    devs <- (syn_mtx_full[, cls, drop=FALSE] - y)/denom
    y <- replicate(lents[x], orig_cnt[, x])
    cnt_full[, cls] <- y + devs*y
  }
  stopifnot(sum(is.na(cnt_full)) == 0) # Validate the assignments
  stopifnot(sum(cnt_full < 0) == 0) # Validate the assignments
  # ===================================
  cnt_full <- as(cnt_full, "dgCMatrix")
  return(cnt_full)
}


# Generates the integrated matrix
ExtendSynthesizedMatrix <- function(syn_mtx, genes, invar_mtx, invarftrs, synt_cells_nm){
  # ===================================
  data_full <- matrix(NA, nrow=length(genes), ncol=ncol(syn_mtx))
  rownames(data_full) <- genes
  colnames(data_full) <- colnames(syn_mtx)
  # ===================================
  # Assign synthesized variable genes data
  indcs <- match(rownames(syn_mtx), rownames(data_full))
  suppressWarnings({
    data_full[indcs, ] <- syn_mtx
  })
  # ===================================
  # Initialize progress bar
  nl <- length(synt_cells_nm)
  pb <- initialize_bar(totbar=nl, wdth=66)
  # ===================================
  # Add invariant genes data to the synthetic matrix, ensuring correct cell name alignment
  indcs <- match(invarftrs, rownames(data_full))
  lents <- lengths(synt_cells_nm)
  for (x in names(synt_cells_nm)) {
    pb$tick() # Update progress bar
    data_full[indcs, synt_cells_nm[[x]]] <- replicate(lents[x], invar_mtx[, x])
  }
  stopifnot(sum(is.na(data_full)) == 0) # Validate the assignments
  # ===================================
  return(data_full)
}


# Generates a list of synthesized cell barcodes and their original counterparts
FetchSynthesizedCells <- function(synt_cells) {
  return(split(synt_cells, gsub("_synth.*$", "", synt_cells)))
}


# Extend the original cells metadata to the integrated object
ExtendMetaData <- function(mtd, mtx){
  mtd[] <- lapply(mtd, function(x) if (is.factor(x)) as.character(x) else x)
  # ===================================
  # Subset the original dataframe based on non-synthetic column names
  mtd_synt <- mtd[sub("_synth.*$", "", colnames(mtx)), ]
  rownames(mtd_synt) <- colnames(mtx)
  return(mtd_synt)
}


# Function to convert list to matrix with unique column names
convert_list_to_matrix <- function(cell_list) {
  # Extract the cell names and make them unique
  cell_names <- paste0(names(cell_list), "_synth")
  unique_cell_names <- make.unique(cell_names)
  # Combine the list into a matrix, with each list element becoming a column
  mat <- do.call(cbind, lapply(cell_list, as.numeric))
  # Assign the unique names to the matrix columns
  colnames(mat) <- unique_cell_names
  return(mat)
}



# calculate the deviation from the originals
DevScGFT <- function(object) {
  # ===================================
  suppressWarnings({
    cells_synt <- colnames(object)[stringr::str_which(colnames(object), paste0("_synth"))]
    syn_mtx <- as.matrix(Seurat::GetAssayData(subset(object, cells=cells_synt), assay="RNA", layer="data"))
  })
  # ===================================
  suppressWarnings({
    cells_orig <- setdiff(colnames(object), cells_synt)
    orig_mtx <- as.matrix(Seurat::GetAssayData(subset(object, cells=cells_orig), assay="RNA", layer="data"))
  })
  # ===================================
  # print(dim(synts))
  # print(dim(origs))
  # ===================================
  # get the data
  stopifnot(all(rownames(orig_mtx) == rownames(syn_mtx)))
  # ===================================
  # get the synthesized cells names
  synt_cells_nm <- FetchSynthesizedCells(colnames(syn_mtx))
  originalCells <- names(synt_cells_nm)
  stopifnot(length(originalCells) == length(unique(originalCells)))
  # ===================================
  # Initialize progress bar
  nl <- length(originalCells)
  pb <- initialize_bar(totbar=nl, wdth=66)
  lents <- lengths(synt_cells_nm)
  # ===================================
  # Loop through each original cell to compute mean deviation
  devs <- sapply(originalCells, function(x){
    pb$tick() # Update progress bar
    # Calculate relative changes for each gene. Compute the row means of changes, ignoring NA values
    denom <- y <- replicate(lents[x], orig_mtx[, x])
    denom[denom == 0] <- 1
    mean(rowMeans((syn_mtx[, synt_cells_nm[[x]], drop=FALSE] - y) / denom))
  })
  # ===================================
  return(devs)
  # ===================================
}

