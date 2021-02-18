sc_compreg <- function(peak.name.intersect.dir,
                       sample1.mat.dir,
                       sample2.mat.dir,
                       motif.match.mat.dir,
                       motif.target.dir,
                       motif.mat.dir,
                       peak.gene.prior.dir,
                       sep.char = '\t',
                       thresh = 0.2,
                       sig.level = 0.05,
                       num.top.tf = 5000,
                       d0.default = 500000,
                       ...) {
    UseMethod("sc_compreg.default")
}
