library(argparser)
library(ConsensusClusterPlus)


argv <- arg_parser('')
argv <- add_argument(argv, "--fpkm",   help="the fpkm file")
argv <- add_argument(argv, '--maxk',   help = "the output directory", default=10)
argv <- add_argument(argv, '--top',    help = "the top ? variant gene", default=1500)
argv <- add_argument(argv, "--outdir", help="the directory of output file")
argv <- parse_args(argv)

fpkmfile <- argv$fpkm
outdir <- argv$outdir
top <- as.numeric(argv$top)
maxk <- as.numeric(argv$maxk)

setwd(outdir)

raw_fpkm<-read.table(fpkmfile, sep="\t", header = T, row.names = 1, check.names=FALSE)
mads <- apply(raw_fpkm,1,mad)

## 选取变异度最大的gene，
## 如果threshold < 1 , 也就是选取前百分之几
## 否则threshold 就是选取前多少个gene

if (top < 1) {
  ## 取变异最大的前20%gene
  n_remain <- round(length(mads) * top)
  mad_fpkm <- raw_fpkm[rev(order(mads))[1:n_remain],]
}else{
  mad_fpkm <- raw_fpkm[rev(order(mads))[1:top],]
}

## 采用中位数中心化
norm_fpkm <- sweep(mad_fpkm, 1, apply(mad_fpkm,1,median,na.rm=T), FUN="-")
norm_fpkm <- as.matrix(norm_fpkm)

## cluster with hc and pearson
results = ConsensusClusterPlus(norm_fpkm, maxK=maxk, reps=1000, pItem=0.9, pFeature=1, title="ConsensusClusterHC", 
	clusterAlg="hc", distance="pearson", seed=1262118388.71279, plot="png", writeTable=T)

## cluster with km and euclidean
results = ConsensusClusterPlus(norm_fpkm, maxK=maxk, reps=1000, pItem=0.9, pFeature=1, title="ConsensusClusterKM", 
	clusterAlg="km", distance="euclidean", seed=1262118388.71279, plot="png", writeTable=T)



