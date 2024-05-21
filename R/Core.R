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
#'                             for synthesis to capture distinct variation modes
#'                             in each cell population, such as cell types or clusters.
#'                             All groups are uniformly scaled for desired expansion.
#' @param    cells             An optional argument specifying the subset of cells
#'                             to be used for synthesis.
#'
#' @return Returns a combined \code{object} of original and synthetic cells. The
#'   synthesized cells are integrated into the \code{data} layer, while their
#'   reversed log-transformed counts are integrated into the \code{counts} layer.
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
                     groups, cells=NULL) {
  # =======================================
  # Check if Seurat package is installed
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running scGFT on a Seurat object requires Seurat.")
  }
  # Check if Seurat object contains RNA assay
  if (!"RNA" %in% Seurat::Assays(object)) {
    stop("Seurat object does not contain RNA assay.")
  }
  # Check if Seurat object contains RNA assay with layers for counts and data
  if (!all(c("counts", "data") %in% names(object$RNA@layers))) {
    stop("Layers for counts or data does not exist. Please run Seurat::NormalizeData() first.")
  }
    # Check if Seurat object contains variable features
  if (object@commands[["NormalizeData.RNA"]]$normalization.method != "LogNormalize") {
    stop("Seurat object should be 'LogNormalize'. Please use `LogNormalize` for 'normalization.method'.")
  }
  scl_fctr <- object@commands[["NormalizeData.RNA"]]$scale.factor
  # Check if Seurat object contains variable features
  if (is.null(Seurat::VariableFeatures(object))) {
    stop("Seurat object does not contain variable features. Please run Seurat::FindVariableFeatures() first.")
  }
  # check if 'groups' exist in the metadata
  if (!groups %in% names(object@meta.data)) {
    stop(paste0("'", groups, "' does not exist in the object metadata."))
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
  # =======================================
  grps <- as.character(object@meta.data[[groups]])
  genes <- rownames(object)
  varftrs <- Seurat::VariableFeatures(object)
  invarftrs <- setdiff(genes, varftrs)
  # =======================================
  if (!is.null(cells)) {
    var_mtx <- orig_dta[varftrs, cells, drop=FALSE]
    invar_mtx <- orig_dta[invarftrs, cells, drop=FALSE]
    grps <- rep(1, length(cells))
  } else {
    var_mtx <- orig_dta[varftrs, ]
    invar_mtx <- orig_dta[invarftrs, ]
  }
  # =======================================
  start_time <- Sys.time()
  syn_mtx <- PerformDIFT(cnt_mtx=var_mtx, groups=grps, nsynth=nsynth, ncpmnts=ncpmnts)
  end_time <- Sys.time()
  dt <- as.numeric(difftime(end_time, start_time, units = "secs"))
  message(paste("Synthesis completed in:", round(dt/60, 2), "min"))
  # ===================================
  message(paste("Integrating data (1/4)"))
  # get the synthesized cells names
  synt_cells_nm <- FetchSynthesizedCells(syn_mtx)
  # ===================================
  message(paste("Integrating data (2/4)"))
  # extend to integrated matrix
  syn_mtx_full <- ExtendSynthesizedMatrix(syn_mtx, genes, invar_mtx, invarftrs, synt_cells_nm)
  # ===================================
  message(paste("Integrating data (3/4)"))
  metadata_synt <- ExtendMetaData(object, syn_mtx_full)
  # ===================================
  message(paste("Integrating data (4/4)"))
  sobj_synt <- SobjMerger(cnt_ls = list(as(orig_dta, "dgCMatrix"), as(syn_mtx_full, "dgCMatrix")),
                          mtd_ls = list(object@meta.data, metadata_synt))
  sobj_synt@assays$RNA$data <- sobj_synt@assays$RNA$counts
  sobj_synt@assays$RNA$counts <- GetCountMatrix(sobj_synt, orig_cnt, orig_dta, genes, synt_cells_nm, syn_mtx_full, scl_fctr)
  Seurat::VariableFeatures(sobj_synt) <- varftrs
  # ===================================
  ids <- grepl("_synth", rownames(sobj_synt@meta.data))
  table(ids)
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
#' @return   synthesis accuracy in percentage.
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
  message(paste("Synthesized cells:", format(results_stats$total_cells, big.mark=",")))
  message(paste("Matching groups:", format(results_stats$matching_clusters, big.mark=",")))
  message(paste("Accuracy (%):", round(results_stats$accuracy_percentage, 2)))
  # ===================================
  message("Calculating deviation from originals...")
  DevScGFT(object)
  # ===================================
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
ModificationGroup <- function(ft_mtx, cmpnts, cellsInGroup) {
  ampltds <- abs(ft_mtx[cmpnts, cellsInGroup, drop = FALSE])
  ampltds <- apply(ampltds, 1, function(x) x/sqrt(sum(x^2)))
  mdfs <- apply(ampltds, 2, function(x) rnorm(1, mean=2, sd=sd(x))) # the previous line flips the dimensions
  return(mdfs)
}


# Sample a modification factor from a Normal Gaussian distribution
ModificationSingle <- function(cmpnts) {
  mdfs <- rnorm(length(cmpnts), mean=2, sd=1)
  return(mdfs)
}


# Performs Discrete and INverse Fourier Transforms
PerformDIFT <- function(cnt_mtx, groups, nsynth, ncpmnts=1) {
  # Extracting components from the input object
  geneNames <- rownames(cnt_mtx)
  cellNames <- colnames(cnt_mtx)
  # originalUMICount <- colSums(cnt_mtx) # Obtain original cell's total UMI count for normalization
  # =======================================
  # Step 1: Calculate the number of cells to synthesize per group
  groupCounts <- table(groups)
  groupProportions <- groupCounts / sum(groupCounts)
  nPerGroup <- round(nsynth * groupProportions)
  # =======================================
  excess_cells <- sum(nPerGroup) - nsynth
  # Adjust the group with the maximum number of cells if there's an excess
  if (excess_cells > 0) {
    max_group <- which.max(nPerGroup)
    nPerGroup[max_group] <- nPerGroup[max_group] - excess_cells
  }
  # Similarly, you can adjust for a deficit in the total (if sum(nPerGroup) < nsynth)
  # by adding missing cells to the group with the least cells
  missing_cells <- nsynth - sum(nPerGroup)
  if (missing_cells > 0) {
    min_group <- which.min(nPerGroup)
    nPerGroup[min_group] <- nPerGroup[min_group] + missing_cells
  }
  # =======================================
  # geneMatrix: Matrix with gene names as rows and cell barcodes as columns
  # groups: Vector containing groups IDs for each cell
  # Step 1: Apply Fourier Transform to each cell
  message(paste("Discrete fourier transform..."))
  ft_mtx <- apply(cnt_mtx, 2, function(cell) fft(cell))
  rownames(ft_mtx) <- paste0("cmpnt_",1:nrow(cnt_mtx))
  len_ft <- nrow(ft_mtx)
  # =======================================
  message(paste("Inverse fourier transform..."))
  message(paste("synthesizing", format(sum(nPerGroup), big.mark=","), "cells..."))
  # Step 2: Generate synthetic cells for each group
  # set.seed(12345)
  synthMatrix_ls <- list()
  cell_cnt <- 0
  for (groupID in names(sort(nPerGroup, decreasing=T))) {
    # get cells in the group
    cellsInGroup <- cellNames[groups == groupID]
    # decide modification function
    if (length(cellsInGroup) > 1) {
      ModificationFunc <- function(cmpnts) ModificationGroup(ft_mtx, cmpnts, cellsInGroup)
    } else {
      ModificationFunc <- function(cmpnts) ModificationSingle(cmpnts)
    }
    # Randomly select a cell from the group
    n_synt <- nPerGroup[names(nPerGroup) == groupID]
    selectedCells <- sample(cellsInGroup, n_synt, replace=TRUE)
    synthMatrix <- lapply(selectedCells, function(x){
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
      # syn_cell <- syn_cell * (originalUMICount[x] / sum(syn_cell)) # Scale UMI counts
      return(syn_cell)
    })
    names(synthMatrix) <- selectedCells
    synthMatrix_ls <- c(synthMatrix_ls, synthMatrix)
    cell_cnt <- cell_cnt + length(synthMatrix)
    message(paste(format(cell_cnt, big.mark=","), "cells synthesized..."))
  }
  stopifnot(length(synthMatrix_ls) == nsynth)
  stopifnot(unique(lengths(synthMatrix_ls)) == length(geneNames))
  synthMatrix <- convert_list_to_matrix(synthMatrix_ls)
  rownames(synthMatrix) <- geneNames
  synthMatrix[synthMatrix < 0] <- 0 # Apply ReLU
  syn_mtx <- as(synthMatrix, "dgCMatrix")
  # =======================================
  return(syn_mtx)
}


# Generated the raw count matrix by reversing log-transformed data
GetCountMatrix <- function(sobj_synt, orig_cnt, orig_dta, genes, synt_cells_nm, syn_mtx_full, scl_fctr) {
  # ===================================
  cnt_full <- matrix(NA, nrow=length(genes), ncol=ncol(sobj_synt))
  rownames(cnt_full) <- genes
  colnames(cnt_full) <- colnames(sobj_synt)
  # ===================================
  indcs <- match(colnames(orig_cnt), colnames(cnt_full))
  suppressWarnings({
    cnt_full[, indcs] <- orig_cnt
  })
  # ===================================
  # Initialize progress bar
  nl <- length(synt_cells_nm)
  pb <- initialize_bar(totbar=nl, wdth=66)
  # ===================================
  # Add count matrix to the full synth matrix, ensuring correct cell name alignment
  # orig_tot_umi <- colSums(orig_cnt)/scl_fctr # Obtain original cell's total UMI count for normalization
  # for (x in names(synt_cells_nm)) {
  #   pb$tick() # Update progress bar
  #   y <- synt_cells_nm[[x]]
  #   cnt_full[, y] <- expm1(syn_mtx_full[, y]) * orig_tot_umi[x]
  # }
  for (x in names(synt_cells_nm)) {
    pb$tick() # Update progress bar
    cls <- synt_cells_nm[[x]]
    devs <- apply(syn_mtx_full[, cls, drop=FALSE], 2, function(y) RelativeChange(orig_dta[, x, drop=FALSE], y))
    cnt_full[, cls] <- apply(devs, 2, function(y) revRelativeChange(orig_cnt[, x, drop=FALSE], y))
  }
  stopifnot(sum(is.na(cnt_full)) == 0) # Validate the assignments
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
    data_full[indcs, ] <- as.matrix(syn_mtx)
  })
  # ===================================
  # Initialize progress bar
  nl <- length(synt_cells_nm)
  pb <- initialize_bar(totbar=nl, wdth=66)
  # ===================================
  # Add invariant genes data to the synthetic matrix, ensuring correct cell name alignment
  indcs <- match(invarftrs, rownames(data_full))
  for (x in names(synt_cells_nm)) {
    pb$tick() # Update progress bar
    y <- synt_cells_nm[[x]]
    data_full[indcs, y] <- replicate(length(y), invar_mtx[, x])
  }
  stopifnot(sum(is.na(data_full)) == 0) # Validate the assignments
  # ===================================
  return(data_full)
}


# Generates a list of synthesized cell barcodes and their original counterparts
FetchSynthesizedCells <- function(mtx) {
  syntCells <- colnames(mtx)
  originalCells <- unique(sub("_synth.*$", "", syntCells))
  # ===================================
  # Initialize progress bar
  nl <- length(originalCells)
  pb <- initialize_bar(totbar=nl, wdth=66)
  # ===================================
  synt_cells_nm <- vector(mode="list", length=nl)
  for (i in 1:nl) {
    pb$tick() # Update progress bar
    synt_cells_nm[[i]] <- syntCells[str_which(syntCells, paste0("^",originalCells[i],"_synth"))]
  }
  names(synt_cells_nm) <- originalCells
  return(synt_cells_nm)
}


# Extend the original cells metadata to the integrated object
ExtendMetaData <- function(sobj, mtx){
  mtd <- sobj@meta.data
  mtd[] <- lapply(mtd, function(x) if (is.factor(x)) as.character(x) else x)
  col_names <- colnames(mtd)
  cols <- vector("list", length = length(col_names))
  names(cols) <- col_names
  # ===================================
  for (i in seq_along(cols)) {
    cols[[i]] <- vector(length=ncol(mtx))
  }
  # ===================================
  metadata_synt <- as.data.frame(cols)
  rownames(metadata_synt) <- colnames(mtx)
  # ===================================
  # Initialize progress bar
  nl <- ncol(mtx)
  pb <- initialize_bar(totbar=nl, wdth=66)
  # ===================================
  for (x in colnames(mtx)) {
    pb$tick() # Update progress bar
    metadata_synt[x, ] <- mtd[sub("_synth.*$", "", x), ]
  }
  return(metadata_synt)
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
  synt_cells_nm <- FetchSynthesizedCells(syn_mtx)
  originalCells <- names(synt_cells_nm)
  stopifnot(length(originalCells) == length(unique(originalCells)))
  # ===================================
  # Initialize progress bar
  nl <- length(originalCells)
  pb <- initialize_bar(totbar=nl, wdth=66)
  # ===================================
  # Loop through each original cell to compute mean deviation
  dev_mtx <- sapply(originalCells, function(x){
    pb$tick() # Update progress bar
    # Calculate relative changes for each gene. Compute the row means of changes, ignoring NA values
    mean(rowMeans(apply(syn_mtx[, synt_cells_nm[[x]], drop=FALSE], 2, function(y) RelativeChange(orig_mtx[, x, drop=FALSE], y)), na.rm=TRUE))
  })
  # ===================================
  message(paste("Deviation (%):", round(mean(dev_mtx)*100, 2), "+/-", round(sd(dev_mtx)*100, 2)))
  # ===================================
}


# Function to safely calculate relative changes to avoid division by zero
RelativeChange <- function(a, b) {
  change <- (b - a) / ifelse(a == 0, 1, a)
  return(change)
}

# Function to safely calculate relative changes to avoid division by zero
revRelativeChange <- function(a, b) {
  return(a+b*a)
}
