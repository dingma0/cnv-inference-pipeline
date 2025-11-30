#==== Parse arguments ====
#packageVersion('numbat') ‘1.5.0’
library(optparse)
options(stringsAsFactors = F)
option_list = list(make_option("--countmat", default = NULL),
                   make_option("--alleledf", default = NULL),
                   make_option("--ncores", default = NULL),
                   make_option("--outdir", default = NULL)
                   )
args = parse_args(OptionParser(option_list = option_list))

library(numbat)
library(data.table)

obj = readRDS(args$countmat)
count_mat = obj@assays$RNA$counts
df_allele = fread(args$alleledf)

out = run_numbat(
    count_mat,
    ref_hca,
    df_allele,
    genome = "hg38",
    t = 1e-5,
    ncores = as.integer(args$ncores),
    plot = TRUE,
    out_dir = args$outdir,
    check_convergence = TRUE,
    multi_allelic = FALSE,
    skip_nj = TRUE,
    max_iter = 1
)
