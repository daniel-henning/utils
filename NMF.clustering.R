suppressMessages(if(!require(argparser)){BiocManager::install('argparser')})
suppressMessages(if(!require(NMF)){BiocManager::install('NMF')})
suppressMessages(if(!require(IntNMF)){BiocManager::install('IntNMF')})
## parse arguments
argv <- arg_parser('NMF cluster')
argv <- add_argument(argv, '--fpkm', help = "all the sample fpkm file")
argv <- add_argument(argv, '--outdir', help = "the output directory")
argv <- add_argument(argv, '--maxk', help = "the output directory", default=10)
argv <- add_argument(argv, '--top', help = "the top ? variant gene", default=1500)
argv <- parse_args(argv)

## get file
fpkmfile <- argv$fpkm
outdir <- argv$outdir
top <- as.numeric(argv$top)
maxk <- as.numeric(argv$maxk)

setwd(outdir)

raw_fpkm<-read.table(fpkmfile, sep="\t", header = T, row.names = 1, check.names=FALSE)

mads <- apply(raw_fpkm,1,mad)

if (top < 1) {
  ## select top % genes
  n_remain <- round(length(mads) * top)
  mad_fpkm <- raw_fpkm[rev(order(mads))[1:n_remain],]
}else{
  mad_fpkm <- raw_fpkm[rev(order(mads))[1:top],]
}

# estimate rank k value and plot the results
fit <- nmf(mad_fpkm, 2:maxk, nrun=10, seed=123456)
pdf(paste0(outdir, "/NMF_rank_survey.pdf"), onefile = FALSE)
plot(fit)
dev.off()

# runing clustering analysis by differential rank k.
for (i in 2:maxk) {
  nmf.each <- nmf(mad_fpkm, i, nrun=10, seed=123456)
  pdf(paste0(outdir, "/consensusNMF_heatmap_k", i, ".pdf"), onefile = FALSE)
  consensusmap(nmf.each, tracks = "consensus")
  dev.off()
  
  outfile = paste0(outdir, "/consensusNMF_matrix_k", i, ".xls")
  mat <- nmf.each@consensus
  res <-  cbind(Sample=row.names(mat), mat)
  write.table(res, file = outfile, sep= "\t", row.names = FALSE, quote = FALSE)
  
  cluster = paste0(outdir, "/consensusNMF_group_k", i, ".xls")
  clu <- predict(nmf.each)
  write.table(clu, file=cluster, sep="\t", quote=FALSE, col.names=FALSE)
  
}

##### chosse the optimal K.

dat <- t(mad_fpkm)

nmf.opt.k(dat = dat, n.runs = 30, n.fold = 5, k.range = 2:maxk, result = TRUE, make.plot = TRUE, progress = TRUE, st.count = 10, maxiter = 100)

file.rename('Rplots.pdf', 'opt.k.pdf')


