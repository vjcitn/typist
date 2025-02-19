#' acquire the python object corresponding to a celltypist model
#' @param mname character(1) model name with .pkl suffix, or path to a pickle file holding a model
#' @note This hands back a live reference to a python object, which does
#' not respect the typical basilisk protocol.  A safer version would extract
#' atomic R data and package it as the return value.
#' @examples
#' pa = system.file(file.path("pickle", "Human_PF_Lung.pkl"), package="typistR")
#' model_desc(pa)
#' @export
model_desc = function (mname="Human_Lung_Atlas.pkl") 
{
    proc <- basiliskStart(typistR:::bsklenv)
    on.exit(basiliskStop(proc))
    basiliskRun(proc, function(mname) {
        ct = import("celltypist")
        chk = ct$models$models_description()
        if (!file.exists(mname))
            if (!(mname %in% chk$model)) stop("mname not in model collection")
        ct$models$Model$load(model = mname)
    }, mname=mname)
}
