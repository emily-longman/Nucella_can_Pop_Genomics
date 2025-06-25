#### Moments admixture with Nucella -- for trios
# module load Rgeospatial
library(genomalicious)
library(SeqArray)
library(foreach)
library(sp)
library(poolfstat)


args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
pairsF=args[2]

pairs = fread(pairsF)
#1.Create output the directory
system(paste("mkdir", "/data/processed/pop_structure/dadi_objects_trios", sep = " "))
system(paste("mkdir", "/data/processed/pop_structure/L_meta_objects_trios", sep = " "))

#2. Load the pairs data
SAMPLE_INPUT=c(pairs$Parent1[i],pairs$Parent2[i],pairs$Derived[i])
#

#### Read in GDS file
poolobj <- get(load("/data/processed/pop_structure/pooldata.RData"))

sample.names <- poolobj@poolnames

trio.set <-
  pooldata.subset(
    poolobj,
    pool.index = which(sample.names %in% SAMPLE_INPUT),
    min.cov.per.pool = 10,
    max.cov.per.pool = 150,
    min.maf = 0.01,
    return.snp.idx = TRUE,
    verbose = TRUE
  )

### get allele frequency data
ad <- trio.set@refallele.readcount
dp <- trio.set@readcoverage  
dat <- ad/dp

### parse
dim(dat)
colnames(dat) <- trio.set@poolnames

snp.dt = trio.set@snp.info
snp.dt %>% mutate(snp_id = rownames(.)) -> snp.dt
####
### convert to format for genomalicious
tf <- TRUE
dat <- as.data.table(dat)
dat[,locus:=snp.dt$snp_id[tf]]
dat[,ref:=snp.dt$RefAllele[tf]]
dat[,alt:=snp.dt$AltAllele[tf]]

datl <- melt(dat, id.vars=c("locus", "ref", "alt"))
setnames(datl, c("variable", "value"), c("POOL", "FREQ"))
######

### Calculate NEff
dp = data.frame(dp)
names(dp) = trio.set@poolnames
colMeans(dp, na.rm = T) %>%
  data.frame(DP=.) %>%
  mutate( sampleId=rownames(.)) ->
  dp.df
####


####
dp.df %>%
  mutate(ninds = 20) %>%
  group_by(sampleId) %>%
  mutate(ne = floor((ninds*DP)/(ninds+DP-1)) )->neff

datl <- merge(datl, neff[,c("sampleId", "ne")], by.x="POOL", by.y="sampleId")
datl[,FREQ:=1-FREQ]
datl <- na.omit(datl)

datl %>%
  filter(POOL == pairs$Parent1[i]) %>%
  mutate(POOL = "A_parent1") ->
  datl1

datl %>%
  filter(POOL == pairs$Parent2[i]) %>%
  mutate(POOL = "B_parent2") ->
  datl2

datl %>%
  filter(POOL == pairs$Derived[i]) %>%
  mutate(POOL = "C_derived") ->
  datl3


rbind(datl1, datl2, datl3) -> datl.relab

sfs_method="probs"

dadi <- dadi_inputs(
  datl.relab,
  type = "freqs",
  popCol = "POOL",
  locusCol = "locus",
  refCol = "ref",
  altCol = "alt",
  freqCol = "FREQ",
  indsCol = "ne",
  freqMethod = sfs_method
)
dadi <- na.omit(dadi) 

L=dim(dadi)[1]
####
### Write dadi file

fn <- paste("dadi_objects_trios/",
            sfs_method, ".",
            pairs[i,]$Parent1, ".",
            pairs[i,]$Parent2, ".",
            pairs[i,]$Derived, ".",
            "delim", sep="")

write.table(dadi,
            file=fn,
            sep="\t", quote=F, row.names=F)

### Metadata
neff %>%
  filter(sampleId == pairs[i,]$Parent1) %>% .$ne -> an1_ne
neff %>%
  filter(sampleId == pairs[i,]$Parent2) %>% .$ne -> an2_ne
neff %>%
  filter(sampleId == pairs[i,]$Derived) %>% .$ne -> Der_ne

fn2 <- paste("L_meta_objects_trios/",
            sfs_method, ".",
            pairs[i,]$Parent1, ".",
            pairs[i,]$Parent2, ".",
            pairs[i,]$Derived, ".",
            "meta", sep="")


write.table(data.frame(an1_ne, an2_ne, Der_ne, L),
            file=fn2,
            sep="\t", quote=F, row.names=F)

