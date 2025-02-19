
bsklenv <- basilisk::BasiliskEnvironment(
    pkgname="typistR", envname="typist_env", packages="pandas==2.0.0",
    pip=c("anndata==0.11.1", "celltypist==1.6.3"))

#' Simple demonstration, not basilisk-compliant, but hands back python reference
#' @importFrom reticulate import
#' @importFrom basilisk BasiliskEnvironment basiliskStart basiliskStop basiliskRun
#' @export
get_ct <- function() {
    proc <- basiliskStart(bsklenv)
    on.exit(basiliskStop(proc))
    basiliskRun(proc, function() {
        # read in 'SpatialData' from .zarr store
        import("celltypist")
        })}

