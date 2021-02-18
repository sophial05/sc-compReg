sc_compreg.default <- function(peak.name.intersect.dir,
                       sample1.mat.dir,
                       sample2.mat.dir,
                       motif.mat.dir,
                       motif.target.dir,
                       peak.gene.prior.dir,
                       sep.char = '\t',
                       thresh = 0.2,
                       sig.level = 0.05,
                       num.top.tf = 5000,
                       d0.default = 500000) {
    pni.file <- comp_reg_preprocess(peak.name.intersect.dir, token=sep.char)
    s2 <- readMat(sample2.mat.dir)
    s1 <- readMat(sample1.mat.dir)
    clust.profile.output <- cluster_profile(s1$O1,
                                            s1$E1,
                                            s1$O1.idx,
                                            s1$E1.idx,
                                            s1$Symbol1,
                                            s1$PeakName1,
                                            s2$O2,
                                            s2$E2,
                                            s2$O2.idx,
                                            s2$E2.idx,
                                            s2$Symbol2,
                                            s2$PeakName2,
                                            pni.file$vo,
                                            pni.file$vt)

    subpop.link.output <- subpopulation_link(clust.profile.output$E1.mean,
                                             clust.profile.output$E2.mean,
                                             clust.profile.output$O1.mean,
                                             clust.profile.output$O2.mean)

    mfbs.output <- mfbs(clust.profile.output$elem.name,
                         motif.target.dir,
                         motif.mat.dir)

    compreg.output <- compreg(clust.profile.output$symbol,
                              mfbs.output$tf.binding,
                              clust.profile.output$elem.name,
                              peak.gene.prior.dir,
                              s1$E1,
                              s1$E1.idx,
                              s2$E2,
                              s2$E2.idx,
                              clust.profile.output$O1.mean,
                              clust.profile.output$O2.mean,
                              subpop.link.output$match,
                              thresh,
                              sig.level,
                              num.top.tf,
                              d0.default)

    compreg.output$call <- match.call()
    return(compreg.output)
}
