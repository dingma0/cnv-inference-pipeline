#!/usr/bin/env Rscript

library(optparse)

parser = OptionParser(description='Run SNP pileup and phasing pre-phased SNPs')
parser = add_option(parser, '--label', default = 'subject', type = "character", help = "Individual label. One per run.")
parser = add_option(parser, '--samples', default = 'sample', type = "character", help = "Sample name(s); comma delimited if multiple. All samples must belong to the same individual.")
parser = add_option(parser, '--bams', default = NULL, type = "character",  help = "BAM file(s); one per sample, comma delimited if multiple.")
parser = add_option(parser, '--barcodes', default = NULL, type = "character", help = "Cell barcode file(s); one per sample, comma delimited if multiple.")
parser = add_option(parser, '--gmap', default = NULL, type = "character", help = "Path to genetic map provided by Eagle2 (e.g. Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz)")
parser = add_option(parser, '--snpvcf', default = NULL, type = "character", help = "SNP VCF for pileup")
parser = add_option(parser, '--outdir', default = NULL, type = "character", help = "Output directory")
parser = add_option(parser, '--ncores', default = 1, type = "integer", help = "Number of cores")
parser = add_option(parser, '--UMItag', default = "Auto", type = "character", help = "UMI tag(s) in bam; should be 'Auto' or 'UB' for 10x scRNA-seq, 'None' for scATAC-seq, and XM for Slide-seq; comma delimited if multiple data modalities in different samples.")
parser = add_option(parser, '--cellTAG', default = "CB", type = "character", help = "Cell tag in bam; should be CB for 10x and XC for Slide-seq")
parser = add_option(parser, '--smartseq', default = FALSE, action = 'store_true', help = "Running with SMART-seq mode; Supply a txt file containing directories of BAM files to --bams and a txt file containing cell names to --barcodes (each entry on its own line for both; ordering must match).")
parser = add_option(parser, '--bulk', default = FALSE, action = 'store_true', help = "Running with bulk RNA-seq mode; supply --bams and --samples but not --barcodes.")
args <- parse_args(parser)

suppressPackageStartupMessages({
    library(glue)
    library(stringr)
    library(data.table)
    library(dplyr)
    library(vcfR)
    library(Matrix)
    library(numbat)
})

for (arg in names(args)) {
  message(paste0(arg, ": ", args[[arg]]))
}

if (any(sapply(list(args$bams, args$snpvcf, args$outdir, args$gmap), is.null))) {
    stop('Missing one or more always required arguments: --bams, --snpvcf, --outdir, --gmap')
}

if (args$smartseq) {
    if (any(sapply(list(args$barcodes), is.null))) {
        stop('Missing one or more required arguments for smartseq mode: --barcodes')
    }
    mode = 'SMART-Seq'
} else if (args$bulk) {
    if (any(sapply(list(args$samples), is.null))) {
        stop('Missing one or more required arguments for bulk mode: --samples')
    }
    mode = 'Bulk'
} else if (length(unique(str_split(args$UMItag, ',')[[1]])) > 1) {
    mode = 'Mixed 10X'
} else {
    if (any(sapply(list(args$samples, args$barcodes), is.null))) {
        stop('Missing one or more required arguments for 10x mode: --samples, --barcodes')
    }
    mode = '10X'
}

label = args$label
samples = str_split(args$samples, ',')[[1]]
outdir = args$outdir
bams = str_split(args$bams, ',')[[1]]
if (!is.null(args$barcodes)) {
    barcodes = str_split(args$barcodes, ',')[[1]]
}
n_samples = length(samples)
ncores = args$ncores
gmap = args$gmap
snpvcf = args$snpvcf

# parse UMI tags
if (str_detect(args$UMItag, ',')) {
    UMItag = str_split(args$UMItag, ',')[[1]]
} else {
    UMItag = rep(args$UMItag, n_samples)
}

# check that every UMI tag is valid
valid_UMI_tags = c('Auto', 'UB', 'None', 'XM')
if (!all(UMItag %in% valid_UMI_tags)) {
    stop(glue('Invalid UMI tag(s) provided: {setdiff(UMItag, valid_UMI_tags)}. Valid options are: {paste(valid_UMI_tags, collapse = ", ")}'))
}

cellTAG = args$cellTAG
smartseq = args$smartseq
bulk = args$bulk
genome = ifelse(str_detect(args$gmap, 'hg19'), 'hg19', 'hg38')

message(paste0('Running in ', mode, ' mode'))
message(paste0('Using genome version: ', genome))

## check if files exist
for (bam in bams) {
    if (!file.exists(bam)) {
        stop('BAM file not found')
    }
}

if (!file.exists(snpvcf)) {
    stop('SNP VCF not found')
}


if (!is.null(barcodes)) {
    for (barcode in barcodes) {
        if (!file.exists(barcode)) {
            stop('Barcode file not found')
        }
    }
}

if (!file.exists(gmap)) {
    stop('Genetic map not found')
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(glue('{outdir}/pileup'), showWarnings = FALSE)
for (sample in samples) {
    dir.create(glue('{outdir}/pileup/{sample}'), showWarnings = FALSE)
}

## pileup

cmds = c()

if (bulk) {

    for (i in 1:n_samples) {

        bam_file = glue('{outdir}/pileup/{sample}/bam_path.tsv')
        sample_file = glue('{outdir}/pileup/{sample}/sample.tsv')

        fwrite(list(bams[i]), bam_file)
        fwrite(list(samples[i]), sample_file)
        
        cmd = glue(
            'cellsnp-lite', 
            '-S {bam_file}',
            '-i {sample_file}',
            '-O {outdir}/pileup/{samples[i]}',
            '-R {snpvcf}', 
            '-p {ncores}',
            '--minMAF 0',
            '--minCOUNT 2',
            '--UMItag None',
            '--cellTAG None',
            .sep = ' ')

        cmds = c(cmds, cmd)

    }

} else if (smartseq) {

    cmd = glue(
            'cellsnp-lite', 
            '-S {bams}',
            '-i {barcodes}',
            '-O {outdir}/pileup/{samples}',
            '-R {snpvcf}', 
            '-p {ncores}',
            '--minMAF 0',
            '--minCOUNT 2',
            '--UMItag None',
            '--cellTAG None',
            .sep = ' ')

    cmds = c(cmd)

} else {
    
    for (i in 1:n_samples) {
        
        cmd = glue(
            'cellsnp-lite', 
            '-s {bams[i]}',
            '-b {barcodes[i]}',
            '-O {outdir}/pileup/{samples[i]}',
            '-R {snpvcf}', 
            '-p {ncores}',
            '--minMAF 0',
            '--minCOUNT 2',
            '--UMItag {UMItag[i]}',
            '--cellTAG {cellTAG}',
            .sep = ' ')

        cmds = c(cmds, cmd)

    }

}

cat('Running pileup\n')

script = glue('{outdir}/run_pileup.sh')

list(cmds) %>% fwrite(script, sep = '\n')

system(glue('chmod +x {script}'))

system(glue('sh {script} 2>&1 | tee {outdir}/pileup.log'), intern = FALSE)

## Generate allele count dataframe
cat('Generating allele count dataframes\n')

if (genome == 'hg19') {
    gtf = gtf_hg19
} else {
    gtf = gtf_hg38
}

genetic_map = fread(gmap) %>% 
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

for (sample in samples) {
    
    # read in phased VCF
    vcf_phased <- fread(snpvcf, skip = '#CHROM') %>%
        dplyr::rename(CHROM = `#CHROM`) %>%
        dplyr::mutate(CHROM = str_remove(CHROM, 'chr')) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM)))

    last_col <- names(vcf_phased)[ncol(vcf_phased)]
    vcf_phased <- mutate(vcf_phased, !!label := stringr::str_remove(get(last_col), ":.*"))

    pu_dir = glue('{outdir}/pileup/{sample}')

    # pileup VCF
    vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM') %>% 
        rename(CHROM = `#CHROM`) %>%
        mutate(CHROM = str_remove(CHROM, 'chr'))

    # count matrices
    AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
    DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))

    cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)

    df = numbat:::preprocess_allele(
        sample = label,
        vcf_pu = vcf_pu,
        vcf_phased = vcf_phased,
        AD = AD,
        DP = DP,
        barcodes = cell_barcodes,
        gtf = gtf,
        gmap = genetic_map
    )
    
    fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
}

cat('All done!\n')