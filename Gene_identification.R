#library(package=ISLR, lib.loc = "~/Documents/probStats/finals/")
library(package=ISLR, lib.loc = "G:/")
# this function outputs an R value, below which all the Hypothesis are rejected
# FDR <- function(pvalues, q) {
#   m <- length(pvalues)
#   R <- NULL
#   pvalues <- sort.int(pvalues, decreasing = FALSE)
#   
#   for (i in 1:m) {
#     if (pvalues[i] < (i*q/m)) {
#       R <- i
#     }
#   }
#   list(R = R)
# }

FDR <- function(pvalues) {
  alpha <- seq(0, 0.2, 0.001)
  alpha <- sort(alpha, decreasing = TRUE)
  R <- 101
  x <- 0
  while(R>100) {
    x <- x+1
    a <- alpha[x]
    ps <- sort(pvalues, decreasing = FALSE)
    m <- length(pvalues)
    l <- c((1:m)*a/m)
    com <- ps - l
    R <- which(com<=0)[length(which(com<=0))]
    if(length(R) == 0) {
      R <- 0
    }
  }
  Tvalue <- ps[R]
  rej <- which(pvalues<Tvalue)
  print(rej)
  # return(rej)
  # print(R)
  # print(x)
}

PValue_Generator <- function(lis) {
  mu_list <- apply(lis[2:nrow(lis),], 2, mean)
  sd_list <- apply(lis[2:nrow(lis),], 2, sd)
  tlist <- mu_list/(sd_list/sqrt(nrow(lis) - 1))
  tlist <- na.omit(tlist)
  plist.suppressed <- pt(tlist, nrow(lis)-2)
  plist.active <- 1 - pt(tlist, nrow(lis)-2)
  list(plist.active = plist.active, plist.suppressed = plist.suppressed)
}

Main_function <- function() {
  breast.mat <- matrix(ncol = 6830)
  breast.rows <- which(NCI60$labs %in% "BREAST")
  for (i in breast.rows) {
    breast.mat <- rbind(breast.mat, NCI60$data[i,])
  }
  
  cns.mat <- matrix(ncol = 6830)
  cns.rows <- which(NCI60$labs %in% "CNS")
  for (i in cns.rows) {
    cns.mat <- rbind(cns.mat, NCI60$data[i,])
  }
  
  colon.mat <- matrix(ncol = 6830)
  colon.rows <- which(NCI60$labs %in% "COLON")
  for (i in colon.rows) {
    colon.mat <- rbind(colon.mat, NCI60$data[i,])
  }
  
  ka.repro.mat <- matrix(ncol = 6830)
  ka.repro.rows <- which(NCI60$labs %in% "K562A-repro")
  for (i in ka.repro.rows) {
    ka.repro.mat <- rbind(ka.repro.mat, NCI60$data[i,])
  }
  
  kb.repro.mat <- matrix(ncol = 6830)
  kb.repro.rows <- which(NCI60$labs %in% "K562B-repro")
  for (i in kb.repro.rows) {
    kb.repro.mat <- rbind(kb.repro.mat, NCI60$data[i,])
  }
  
  leukemia.mat <- matrix(ncol = 6830)
  leukemia.rows <- which(NCI60$labs %in% "LEUKEMIA")
  for (i in leukemia.rows) {
    leukemia.mat <- rbind(leukemia.mat, NCI60$data[i,])
  }
  
  ma.repro.mat <- matrix(ncol = 6830)
  ma.repro.rows <- which(NCI60$labs %in% "MCF7A-repro")
  for (i in ma.repro.rows) {
    ma.repro.mat <- rbind(ma.repro.mat, NCI60$data[i,])
  }
  
  md.repro.mat <- matrix(ncol = 6830)
  md.repro.rows <- which(NCI60$labs %in% "MCF7D-repro")
  for (i in md.repro.rows) {
    md.repro.mat <- rbind(md.repro.mat, NCI60$data[i,])
  }
  
  melanoma.mat <- matrix(ncol = 6830)
  melanoma.rows <- which(NCI60$labs %in% "MELANOMA")
  for (i in melanoma.rows) {
    melanoma.mat <- rbind(melanoma.mat, NCI60$data[i,])
  }
  
  nsclc.mat <- matrix(ncol = 6830)
  nsclc.rows <- which(NCI60$labs %in% "NSCLC")
  for (i in nsclc.rows) {
    nsclc.mat <- rbind(nsclc.mat, NCI60$data[i,])
  }
  
  ovarian.mat <- matrix(ncol = 6830)
  ovarian.rows <- which(NCI60$labs %in% "OVARIAN")
  for (i in ovarian.rows) {
    ovarian.mat <- rbind(ovarian.mat, NCI60$data[i,])
  }
  
  prostate.mat <- matrix(ncol = 6830)
  prostate.rows <- which(NCI60$labs %in% "PROSTATE")
  for (i in prostate.rows) {
    prostate.mat <- rbind(prostate.mat, NCI60$data[i,])
  }
  
  renal.mat <- matrix(ncol = 6830)
  renal.rows <- which(NCI60$labs %in% "RENAL")
  for (i in renal.rows) {
    renal.mat <- rbind(renal.mat, NCI60$data[i,])
  }
  
  unknown.mat <- matrix(ncol = 6830)
  unknown.rows <- which(NCI60$labs %in% "UNKNOWN")
  for (i in unknown.rows) {
    unknown.mat <- rbind(unknown.mat, NCI60$data[i,])
  }
  
  breast.supp <<- PValue_Generator(breast.mat)$plist.suppressed
  breast.active <<- PValue_Generator(breast.mat)$plist.active
  print("breast.suppressed")
  FDR(breast.supp)
  print("breast.active")
  FDR(breast.active)
  
  cns.supp <<- PValue_Generator(cns.mat)$plist.suppressed
  cns.active <<- PValue_Generator(cns.mat)$plist.active
  print("cns.suppressed")
  FDR(cns.supp)
  print("cns.active")
  FDR(cns.active)
  
  colon.supp <<- PValue_Generator(colon.mat)$plist.suppressed
  colon.active <<- PValue_Generator(colon.mat)$plist.active
  print("colon.suppressed")
  FDR(colon.supp)
  print("colon.active")
  FDR(colon.active)
  
  leukemia.supp <<- PValue_Generator(leukemia.mat)$plist.suppressed
  leukemia.active <<- PValue_Generator(leukemia.mat)$plist.active
  print("leukemia.suppressed")
  FDR(leukemia.supp)
  print("leukemia.active")
  FDR(leukemia.active)
  
  melanoma.supp <<- PValue_Generator(melanoma.mat)$plist.suppressed
  melanoma.active <<- PValue_Generator(melanoma.mat)$plist.active
  print("melanoma.suppressed")
  FDR(melanoma.supp)
  print("melanoma.active")
  FDR(melanoma.active)
  
  nsclc.supp <<- PValue_Generator(nsclc.mat)$plist.suppressed
  nsclc.active <<- PValue_Generator(nsclc.mat)$plist.active
  print("nsclc.suppressed")
  FDR(nsclc.supp)
  print("nsclc.active")
  FDR(nsclc.active)
  
  ovarian.supp <<- PValue_Generator(ovarian.mat)$plist.suppressed
  ovarian.active <<- PValue_Generator(ovarian.mat)$plist.active
  print("ovarian.suppressed")
  FDR(ovarian.supp)
  print("ovarian.active")
  FDR(ovarian.active)
  
  prostate.supp <<- PValue_Generator(prostate.mat)$plist.suppressed
  prostate.active <<- PValue_Generator(prostate.mat)$plist.active
  print("prostate.suppressed")
  FDR(prostate.supp)
  print("prostate.active")
  FDR(prostate.active)
  
  renal.supp <<- PValue_Generator(renal.mat)$plist.suppressed
  renal.active <<- PValue_Generator(renal.mat)$plist.active
  print("renal.suppressed")
  FDR(renal.supp)
  print("renal.active")
  FDR(renal.active)
  
}