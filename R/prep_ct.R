#' list available celltypist models
#' @import basilisk
#' @export
ct_models = function() {
  proc <- basiliskStart(bsklenv)
  on.exit(basiliskStop(proc))
  basiliskRun(proc, function() {
        ct = import("celltypist")
        ct$models$models_description()
        })
  }
 
#' perform predictions
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import sparseMatrixStats
#' @import zellkonverter
#' @param sce SingleCellExperiment instance, assumed to have counts as assay component
#' @param model character(1) name of a celltypist model with .pkl suffix
#' @param majority_voting logical(1), if TRUE, a super-clustering procedure is run
#' and voting across clustering outcomes is used to assign cell type
#' @return A SingleCellExperiment with additional colData columns `predicted_labels` and
#' `conf_score` from celltypist
#' @note The input is converted to AnnData h5ad by zellkonverter, and the output
#' of celltypist prediction is bound to the h5ad, which is then converted back to SingleCellExperiment.
#' @examples
#' z = scRNAseq::ZilionisLungData()
#' kp = c("Ciliated cells", "Club cells", "Endothelial cells", "Fibroblasts")
#' kpi = which(z$`Major cell type` %in% kp)
#' zlim = z[, kpi]
#' r = ct_run(zlim, model="Human_Lung_Atlas.pkl")
#' r
#' # examine details added on major cell types with celltypist
#' tail(sort(table(r$`Major.cell.type`, r$predicted_labels)["Fibroblasts",]))
#' tail(sort(table(r$`Major.cell.type`, r$predicted_labels)["Endothelial cells",]))
#' # examine concordance between Zilionis minor subset and celltypist label
#' tail(sort(table(r$`Minor.subset`, r$predicted_labels)["Endo1_PECAM1/CLEC14A/CD34",]))
#' tail(sort(table(r$`Minor.subset`, r$predicted_labels)[,"EC general capillary"]))
#' @export
ct_run = function(sce, model='Immune_All_Low.pkl', majority_voting=FALSE) {
  message("checking for 0-count cells")
  cs = sparseMatrixStats::colSums2(SummarizedExperiment::assay(sce)) #remove zero count cells
  zz = which(cs==0)
  if (length(zz)>0) sce = sce[,-zz]
  proc <- basiliskStart(bsklenv)
  on.exit(basiliskStop(proc))
  basiliskRun(proc, function(sce, model, majority_voting) {
        message("converting SCE")
        ad = zellkonverter::SCE2AnnData(sce)
        message("done.")
        sc = import("scanpy")
        ct = import("celltypist")
        sc$pp$normalize_total(ad, target_sum=10000)  # required by celltypist
        sc$pp$log1p(ad)
        ctm = ct$models$Model$load(model = model)
        p = ct$annotate(ad, model=ctm, majority_voting=majority_voting)
        adat = p$to_adata()
        zellkonverter::AnnData2SCE(adat)
     }, sce=sce, model=model, majority_voting=majority_voting)
}


#>>> preds = celltypist.annotate(zpl, model=mipf)
#🔬 Input data has 173881 cells and 41861 genes
#🔗 Matching reference genes in the model
#🧬 3814 features used for prediction
#⚖️ Scaling input data
#🖋️ Predicting labels
#✅ Prediction done!
#>>> preds
#CellTypist prediction result for 173881 query cells
#    predicted_labels: data frame with 1 column ('predicted_labels')
#    decision_matrix: data frame with 173881 query cells and 38 cell types
#    probability_matrix: data frame with 173881 query cells and 38 cell types
#    adata: AnnData object referred
#>>> preds['adata']
#Traceback (most recent call last):
#  File "<stdin>", line 1, in <module>
#TypeError: 'AnnotationResult' object is not subscriptable
#>>> preds.adata
#AnnData object with n_obs × n_vars = 173881 × 41861
#    obs: 'Library', 'Barcode', 'Patient', 'Tissue', 'Used', 'Barcoding emulsion', 'Total counts', 'Percent counts from mitochondrial genes', 'Most likely LM22 cell type', 'Major cell type', 'Minor subset', 'used_in_NSCLC_all_cells', 'used_in_NSCLC_and_blood_immune', 'used_in_NSCLC_immune', 'used_in_NSCLC_non_immune', 'used_in_blood'
#    uns: 'X_name', 'log1p'
#    obsm: 'SPRING_NSCLC_all_cells', 'SPRING_NSCLC_and_blood_immune', 'SPRING_NSCLC_immune', 'SPRING_NSCLC_non_immune', 'SPRING_blood'
##>>> 
