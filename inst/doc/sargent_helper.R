library("sargent")
library("readxl")

# Use 'cellstatePipe' function below to annotate cells using sargent.
# sobj should contain a 'count' layer as well as 'RNA_nn' attribute.
sobj <- cellstatePipe(geneset_path, sobj)

# Sargent Helper function
cellstatePipe <- function(geneset_path, sobj) {
  # ===================================
  set.seed(42)
  # ===================================
  # readin genesets
  exls <- excel_sheets(file.path(geneset_path, "genesets.xlsx"))
  # positive genes
  gpos <- read_excel(path=file.path(geneset_path, "genesets.xlsx"), 
                     sheet="positive", col_names=TRUE, col_types="text")
  gpos <- lapply(as.list(gpos), function(x){
    x <- x[!is.na(x)]
  })
  gpos[lengths(gpos) == 0] <- NULL
  names(gpos) <- toupper(names(gpos))
  # negative genes
  if ("negative" %in% exls) {
    gneg <- read_excel(path=file.path(geneset_path, "genesets.xlsx"), 
                       sheet="negative", col_names=TRUE, col_types="text")
    gneg <- lapply(as.list(gneg), function(x){
      x <- x[!is.na(x)]
    })
    gneg[lengths(gneg) == 0] <- NULL
    names(gneg) <- toupper(names(gneg))
  } else { gneg <- NULL }
  # ===================================
  # sargent annotation
  sargent_anot <- sargentAnnotation(gex=GetAssayData(sobj, assay="RNA", layer="counts"), 
                                    gene.sets=gpos, 
                                    gene.sets.neg=gneg,
                                    adjacent.mtx=attr(sobj, which="graphs")[["RNA_nn"]])
  sargent_idents <- fetchAssignment(sargent_anot)
  print(sort(table(sargent_idents)))
  # ===================================
  # add new metadata
  sobj <- AddMetaData(
    object = sobj,
    metadata = sargent_idents[rownames(sobj@meta.data)],
    col.name = 'sargent_celltype'
  ) 
  # ===================================  
  return(sobj)
}