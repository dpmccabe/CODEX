segment <- function(Y_qc, Yhat, optK, K, sampname_qc, ref_qc, chr, lmax, mode) {
  lmax <- lmax - 1
  
  res <- adply(1:ncol(Y_qc), 1, function(sampno) {
    message("Segmenting sample ", sampno, ": ", sampname_qc[sampno], ".")
    y <- Y_qc[, sampno]
    yhat <- Yhat[[which(K == optK)]][, sampno]
    num <- length(y)
    y <- c(y, rep(0, lmax))
    yhat <- c(yhat, rep(0, lmax))
    i <- rep(1:num, rep(lmax, num))
    j <- rep(1:lmax, num) + i
    yact <- rep(0, length(i))
    lambda <- rep(0, length(i))
    
    for (k in 1:num) {
      yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + lmax)])[-1]
      lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + lmax)])[-1]
    }
    
    i <- i[j <= num]
    yact <- yact[j <= num]
    lambda <- lambda[j <= num]
    j <- j[j <= num]
    yact[lambda < 20] <- 20
    lambda[lambda < 20] <- 20
    
    if (mode == "integer") {
      chat <- round(2 * (yact / lambda))
    } else if (mode == "fraction") {
      chat <- 2 * (yact / lambda)
    }
    
    lratio <- (1 - chat / 2) * lambda + log((chat + 1e-04) / 2.0001) * yact
    chat[chat > 5] <- 5
    
    samp_data <- data.table(samp = sampname_qc[sampno], chr = chr, i, j, yact, lambda, chat, lratio)[lratio > 0, ]
    
    if (sum(lratio > 0) > 0) {
      if (sum(lratio > 0) >= 2) {
        samp_data <- samp_data[order(-lratio)]
        s <- 1
        
        while (s <= (nrow(samp_data))) {
          samp_data <- rbind(samp_data[s], samp_data[i > samp_data[s, j] | j < samp_data[s, i]])
          s <- s + 1
        }
      } else if (sum(lratio > 0) != 1) {
        stop("unexpected sum(lratio > 0)")
      }
      
      samp_data <- samp_data[order(-lratio)]
      loglikeij <- cumsum(samp_data$lratio)

      mBIC <- aaply(1:nrow(samp_data), 1, function(s) {
        tau <- sort(unique(c(samp_data[1:s, c(i, j)], 1, num)))
        mbic <- loglikeij[s]
        mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] - tau[1:(length(tau) - 1)]))
        mbic <- mbic + (0.5 - (length(tau) - 2)) * log(num)
        return(mbic)
      })

      bic_ind <- 1:which.max(mBIC)
      samp_data <- samp_data[bic_ind]
      samp_data$mbic <- mBIC[bic_ind]
    } else {
      stop("unexpected sum(lratio > 0)")
    }
    
    return(samp_data)
  }, .parallel = T, .id = NULL)
  
  res <- subset(res, mbic > 0)
  
  ret_cols <- c("sample_name", "chr", "cnv", "st_bp", "ed_bp", "length_kb", "st_exon", "ed_exon", "raw_cov", "norm_cov", "lratio", "mBIC")

  if (nrow(res) > 0) {
    res$st_bp <- start(ref_qc)[res$i]
    res$ed_bp <- end(ref_qc)[res$j]

    res$cnv <- NA
    res[res$chat < 2, "cnv"] <- "del"
    res[res$chat > 2, "cnv"] <- "dup"
    res$cnv <- as.factor(res$cnv)
    
    res$length_kb <- (res$ed_bp - res$st_bp + 1) / 1000
    
    colnames(res) <- c("sample_name", "chr", "st_exon", "ed_exon", "raw_cov", "norm_cov", "copy_no", "lratio", "mBIC", "st_bp", "ed_bp", "cnv", "length_kb")
    
    return(res[, ret_cols])
  } else {
    return(NULL)
  }
}
