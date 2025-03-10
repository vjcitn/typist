% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_ct.R
\name{ct_run}
\alias{ct_run}
\title{perform predictions}
\usage{
ct_run(sce, model = "Immune_All_Low.pkl", majority_voting = FALSE)
}
\arguments{
\item{sce}{SingleCellExperiment instance, assumed to have counts as assay component}

\item{model}{character(1) name of a celltypist model with .pkl suffix}

\item{majority_voting}{logical(1), if TRUE, a super-clustering procedure is run
and voting across clustering outcomes is used to assign cell type}
}
\value{
A SingleCellExperiment with additional colData columns `predicted_labels` and
`conf_score` from celltypist
}
\description{
perform predictions
}
\note{
The input is converted to AnnData h5ad by zellkonverter, and the output
of celltypist prediction is bound to the h5ad, which is then converted back to SingleCellExperiment.
}
\examples{
z = scRNAseq::ZilionisLungData()
kp = c("Ciliated cells", "Club cells", "Endothelial cells", "Fibroblasts")
kpi = which(z$`Major cell type` \%in\% kp)
zlim = z[, kpi]
r = ct_run(zlim, model="Human_Lung_Atlas.pkl")
r
# examine details added on major cell types with celltypist
tail(sort(table(r$`Major.cell.type`, r$predicted_labels)["Fibroblasts",]))
tail(sort(table(r$`Major.cell.type`, r$predicted_labels)["Endothelial cells",]))
# examine concordance between Zilionis minor subset and celltypist label
tail(sort(table(r$`Minor.subset`, r$predicted_labels)["Endo1_PECAM1/CLEC14A/CD34",]))
tail(sort(table(r$`Minor.subset`, r$predicted_labels)[,"EC general capillary"]))
}
