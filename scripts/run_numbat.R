#==== Parse arguments ====
#packageVersion('numbat') ‘1.5.0’
library(optparse)
options(stringsAsFactors = F)
option_list = list(make_option("--countmat", default = NULL),
                   make_option("--alleledf", default = NULL),
                   make_option("--ncores", default = NULL),
                   make_option("--outdir", default = NULL),
                   make_option("--ref_mat", default = NULL),
                   make_option("--ref_ids", default = NULL)
                   )
args = parse_args(OptionParser(option_list = option_list))

library(numbat)
library(data.table)

obj = readRDS(args$countmat)
if (inherits(obj, "Seurat")) {
    count_mat = obj@assays$RNA$counts
  } else if (inherits(obj, "dgCMatrix")) {
    count_mat = obj
  } else {
    stop("--countmat expects a Seurat object or dgCMatrix")
  }
  
df_allele = fread(args$alleledf)

ref <- readRDS(args$ref_mat)
if (inherits(ref, "Seurat")) {
    ref_mat = ref@assays$RNA$counts
  } else if (inherits(ref, "dgCMatrix")) {
    ref_mat = ref
  } else {
    stop("--ref_mat expects a Seurat object or dgCMatrix")
  }

ref_annot <- read.table(args$ref_ids, sep = '\t', header = T)
ref_internal = aggregate_counts(ref_mat, ref_annot)

out = run_numbat(
    count_mat,
    ref_internal,
    df_allele,
    genome = "hg38",
    t = 1e-5,
    ncores = as.integer(args$ncores),
    plot = TRUE,
    out_dir = args$outdir
)
