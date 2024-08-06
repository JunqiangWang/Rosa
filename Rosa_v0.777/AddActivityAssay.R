#' This function adds the protein 'activity' assay to a Seurat object. Assay name
#' @export
AddActivityAssay <- function(seu,
                             activity.matrix,
                             assay.name="Activity"
                             ) {
  assay.activity <- CreateAssayObject(data = as.matrix(activity.matrix))
  seu[[assay.name]] <- assay.activity
  seu@active.assay <- 'Activity'
  return(seu)
}
