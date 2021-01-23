mfbs.default <- function(TF.name,
                         elem.name,
                         motif.target.path,
                         motif.name,
                         motif.weight,
                         match2,
                         motif.mat.path) {

    if (! is(motif.target.path, "character")) {
        stop('motif.target.path must be character class. Pleaes check your input.')
    }

    if (! is(motif.mat.path, "character")) {
        stop('motif.mat.path must be character class. Please check your input.')
    }

    motif.mat <- readMat(motif.mat.path)

    motif.name <- unlist(motif.mat$motifName, use.names=F)
    motif.weight <- as.numeric(unlist(motif.mat$motifWeight, use.names=F))
    match2 <- unlist(motif.mat$Match2, use.names = F)
    m2.half.idx <- length(match2) / 2
    match2 <- list('a' = match2[1: m2.half.idx],
                   'b' = match2[(m2.half.idx + 1):length(match2)])

    this.call <- match.call()

    motif.file <- mfbs.load(motif.target.path)
    f3 <- motif.file$C3
    d1 <- is.element(motif.file$C1, elem.name)
    f1 <- match(motif.file$C1, elem.name, nomatch = 0)

    d2 <- is.element(motif.file$C2, motif.name)
    f2 <- match(motif.file$C2, motif.name, nomatch = 0)

    t1 <- setdiff(seq(1, length(motif.name), 1), unique(f2))
    f2 <- c(f2[(d1 * d2) == 1], t1)
    t1.len <- length(t1)
    f1 <- c(f1[(d1 * d2) == 1], rep(1, t1.len))
    f3 <- c(f3[(d1 * d2) == 1], rep(0, t1.len))
    t1 <- setdiff(seq(1, length(elem.name), 1), unique(f1))
    f1 <- c(f1, t1)
    t1.len <- length(t1)
    f2 <- c(f2, rep(1, t1.len))
    f3 <- c(f3, rep(0, t1.len))


    motif.binding <- sparseMatrix(dims = c(length(motif.name),length(elem.name)),
                                  i = f2,
                                  j = f1,
                                  x = f3)
    motif.weight.len <- length(motif.weight)


    motif.binding <- mult(Matrix(diag(1 / (motif.weight + 0.1),
                                      nrow = motif.weight.len,
                                      ncol = motif.weight.len), sparse=T), motif.binding)
    motif.binding@x <- log(motif.binding@x)

    tf.name <- intersect(symbol, unique(match2$b))
    tf.name.len <- length(tf.name)
    tf.binding <- matrix(0, tf.name.len, length(elem.name))

    mf1 <- match(match2$a, motif.name, nomatch=0)
    mf2 <- match(match2$b, tf.name, nomatch=0)

    for (i in 1:tf.name.len) {
        print(i)
        a <- which(mf2 == i)
        if (length(a) > 1) {
            tf.binding[i, ] <- max(motif.binding[mf1[a], ])
        } else if (length(a) == 1) {
            tf.binding[i, ] <- motif.binding[mf1[a], ]
        }
    }

    tf.binding <- Matrix(tf.binding, sparse = T)

    output <- list('tf.binding' = tf.binding)
    output$call = this.call
    return(output)
}