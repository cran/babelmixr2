# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

convertDataBack <- function(id, time, amt, ii, evid, cmt, cmtDvid, dvidDvid, linNcmt = 0L, linKa = 0L, neq = 0L, replaceEvid = 5L, zeroDose2 = TRUE) {
    .Call(`_babelmixr2_convertDataBack`, id, time, amt, ii, evid, cmt, cmtDvid, dvidDvid, linNcmt, linKa, neq, replaceEvid, zeroDose2)
}

transDv <- function(inDv, inCmt, cmtTrans, lambda, yj, low, high) {
    .Call(`_babelmixr2_transDv`, inDv, inCmt, cmtTrans, lambda, yj, low, high)
}

