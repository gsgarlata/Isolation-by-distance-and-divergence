


calcDiffGabi<-function (freqs, metric, pops = row.names(freqs), loci = unique(as.matrix(as.data.frame(strsplit(names(freqs), 
                                                                                                 split = ".", fixed = TRUE), stringsAsFactors = FALSE))[1, 
                                                                                                 ]), global = FALSE, bootstrap = FALSE, n.bootstraps = 1000, 
          object = NULL) 
{
  
  if (!metric %in% c("Fst", "Gst", "Jost's D", "Rst")) {
    stop("metric must be Fst, Gst, Rst, or Jost's D")
  }
  if (!all(pops %in% row.names(freqs))) {
    stop("pops must all be in row names of freqs.")
  }
  freqs <- freqs[pops, ]
  loci <- loci[loci != "Genomes"]
  if (metric == "Rst" && (is.null(object) || any(is.na(Usatnts(object)[loci])))) {
    stop("gendata object required with Usatnts for Rst metric")
  }
  cpd <- function(freqs, metric, loci, bootstrap, n.boostraps) {
    if ("Genomes" %in% names(freqs)) {
      genomes <- freqs$Genomes
      names(genomes) <- row.names(freqs)
      GbyL <- FALSE
    }else {
      GbyL <- TRUE
    }
    hets <- array(0, dim = c(length(loci), 4), dimnames = list(loci, 
                                                               c("HT", "HS", "HTest", "HSest")))
    for (L in loci) {
      if (!GbyL) {
        thesegenomes <- genomes
      }
      else {
        thesegenomes <- freqs[, paste(L, "Genomes", sep = ".")]
      }
      thesefreqs <- freqs[, grep(paste("^", L, "\\.", sep = ""), 
                                 names(freqs)), drop = FALSE]
      thesefreqs <- thesefreqs[, names(thesefreqs) != paste(L, 
                                                            "Genomes", sep = "."), drop = FALSE]
      hsByPop <- apply(as.matrix(thesefreqs), 1, function(x) 1 - 
                         sum(x^2))
      if (metric == "Fst") {
        avgfreq <- unlist(lapply(thesefreqs, weighted.mean, 
                                 w = thesegenomes))
        hets[L, "HS"] <- weighted.mean(hsByPop, thesegenomes)
      }
      if (metric %in% c("Jost's D", "Gst")) {
        avgfreq <- colMeans(thesefreqs)
        hets[L, "HS"] <- mean(hsByPop)
      }
      if (metric == "Rst") {
        replen <- Usatnts(object)[L]
        alleles <- sapply(strsplit(grep(paste("^", L, 
                                              "\\.", sep = ""), names(freqs), value = TRUE), 
                                   ".", fixed = TRUE), function(x) x[2])
        alleles <- alleles[alleles != "Genomes"]
        alleles[alleles == "null"] <- 0
        alleles <- as.integer(alleles)
        if (0 %in% alleles) {
          nullfreqs <- thesefreqs[, match(0, alleles)]
          for (pop in 1:dim(thesefreqs)[1]) {
            thesefreqs[pop, ] <- thesefreqs[pop, ]/(1 - 
                                                      nullfreqs[pop])
          }
          thesegenomes <- round(thesegenomes * (1 - nullfreqs))
        }
        totgenomes <- sum(thesegenomes)
        avgfreq <- colMeans(thesefreqs)
        SSalleledistS <- numeric(dim(freqs)[1])
        SSalleledistT <- 0
        
        if(length(alleles)<2) next
        
        for (i in 1:(length(alleles) - 1)) {
          for (j in (i + 1):length(alleles)) {
            if (alleles[i] == 0 || alleles[j] == 0) 
              next
            sqdiff <- (abs(alleles[i] - alleles[j])/replen)^2
            nocc <- thesefreqs[, i] * thesefreqs[, j] * 
              thesegenomes^2
            SSalleledistS <- SSalleledistS + sqdiff * 
              nocc
            SSalleledistT <- SSalleledistT + sqdiff * 
              totgenomes^2 * avgfreq[i] * avgfreq[j]
          }
        }
        hets[L, "HS"] <- mean(SSalleledistS/(thesegenomes * 
                                               (thesegenomes - 1)))
        hets[L, "HT"] <- SSalleledistT/(totgenomes * 
                                          (totgenomes - 1))
      }
      if (metric %in% c("Fst", "Gst", "Jost's D")) {
        hets[L, "HT"] <- 1 - sum(avgfreq^2)
      }
      if (metric %in% c("Jost's D", "Gst")) {
        meanGenomes <- 1/mean(1/thesegenomes)
        hets[L, "HSest"] <- hets[L, "HS"] * meanGenomes/(meanGenomes - 
                                                           1)
        hets[L, "HTest"] <- hets[L, "HT"] + hets[L, "HSest"]/(2 * 
                                                                meanGenomes)
      }
    }
    if (!bootstrap) {
      n.bootstraps <- 1
    }
    result <- numeric(n.bootstraps)
    for (b in 1:n.bootstraps) {
      if (bootstrap) {
        thishets <- hets[sample(loci, replace = TRUE), 
        ]
      }
      else {
        thishets <- hets
      }
      if (metric == "Fst") {
        HT <- mean(thishets[, "HT"])
        HS <- mean(thishets[, "HS"])
        result[b] <- (HT - HS)/HT
      }
      if (metric == "Rst") {
        R <- (thishets[, "HT"] - thishets[, "HS"])/thishets[,"HT"]
        zero_idx<-as.numeric(which(thishets[, "HT"]==0))
        if(length(zero_idx)==0){}else{
        R[zero_idx]<-0
        }
        result[b] <- mean(R)
      }
      if (metric == "Gst") {
        G <- (thishets[, "HTest"] - thishets[, "HSest"])/thishets[, 
                                                                  "HTest"]
        result[b] <- mean(G)
      }
      if (metric == "Jost's D") {
        D <- 2 * (thishets[, "HTest"] - thishets[, "HSest"])/(1 - 
                                                                thishets[, "HSest"])
        result[b] <- mean(D)
      }
    }
    return(result)
  }
  if (global) {
    result <- cpd(freqs, metric, loci, bootstrap, n.bootstraps)
  }
  else {
    if (bootstrap) {
      result <- array(list(), dim = c(length(pops), length(pops)), 
                      dimnames = list(pops, pops))
    }
    else {
      result <- matrix(0, nrow = length(pops), ncol = length(pops), 
                       dimnames = list(pops, pops))
    }
    for (m in 1:length(pops)) {
      for (n in m:length(pops)) {
        thisres <- cpd(freqs[unique(c(m, n)), ], metric, 
                       loci, bootstrap, n.bootstraps)
        if (bootstrap) {
          result[[m, n]] <- result[[n, m]] <- thisres
        }
        else {
          result[m, n] <- result[n, m] <- thisres
        }
      }
    }
  }
  return(result)
}
