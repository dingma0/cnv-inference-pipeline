process HAPLOTYPECALLER {
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(ref_dir)
    val(ref_name)
    each chrom

    output:
    tuple val(sample_id), path("snps/${sample_id}/${sample_id}_${chrom}.vcf.gz"), path("snps/${sample_id}/${sample_id}_${chrom}.vcf.gz.tbi")

    script:
    // TODO: do something about the bam name situation
    """
    REF=$ref_dir/$ref_name
    MEM='$task.memory'
    mkdir tmp
    BAM=\$(basename -s .bai $bai)
    mv $bam \$BAM
    mkdir -p snps/$sample_id

    gatk --java-options "-Xmx\${MEM%% *}g -XX:-UsePerfData" \
    HaplotypeCaller \
        --input \$BAM \
        --output snps/${sample_id}/${sample_id}_${chrom}.vcf.gz \
        --reference \$REF \
        --native-pair-hmm-threads $task.cpus \
        --intervals $chrom \
        --tmp-dir tmp
    """
}

process BCFTOOLS_HC {
    publishDir "$workflow.outputDir", mode:'copy'

    input:
    tuple val(sample_id), path(vcfs), path(tbi)

    output:
    tuple val(sample_id), path("phasing/${sample_id}/${sample_id}_het_snps.bcf"), path("phasing/${sample_id}/${sample_id}_het_snps.bcf.csi")

    script:
    """
    MEM='$task.memory'

    mkdir -p phasing/$sample_id
    bcftools concat -a -Ob --threads $task.cpus $vcfs | \
        bcftools sort -Ob -m \${MEM%% *}G | \
        bcftools view --threads $task.cpus -m2 -M2 -v snps -Ob | \
        bcftools view --threads $task.cpus -g het -Ob -W=csi -o phasing/${sample_id}/${sample_id}_het_snps.bcf
    """
}

process EAGLE2 {
    input:
    tuple val(sample_id), path(bcf), path(csi)
    path(ref_panel_dir)
    path(gmap)
    each chrom

    output:
    tuple val(sample_id), path("phasing/${sample_id}/${sample_id}_${chrom}_phased_snps.bcf")

    script:
    """
    mkdir -p phasing/$sample_id

    eagle \
        --geneticMapFile $gmap \
        --vcfRef ${ref_panel_dir}/${chrom}.genotypes.bcf \
        --vcfTarget $bcf \
        --chrom $chrom \
        --numThreads $task.cpus \
        --vcfOutFormat b \
        --outPrefix phasing/${sample_id}/${sample_id}_${chrom}_phased_snps
    """
}

process BCFTOOLS_EAGLE {
    publishDir "$workflow.outputDir", mode:'copy'

    input:
    tuple val(sample_id), path(bcfs)

    output:
    tuple val(sample_id), path("phasing/${sample_id}/${sample_id}_phased_snps.bcf"), path("phasing/${sample_id}/${sample_id}_phased_snps.bcf.csi")

    script:
    """
    MEM='$task.memory'
    mkdir -p phasing/$sample_id
    mkdir -p tmp

    for BCF in $bcfs; do 
        bcftools index -c --threads $task.cpus \$BCF 
    done

    bcftools concat -a -Ob --threads $task.cpus $bcfs | \
        bcftools sort -Ob -m \${MEM%% *}G -W=csi -T tmp -o phasing/${sample_id}/${sample_id}_phased_snps.bcf 
    """
}

process DOWNLOAD_TEST_DATA {
    publishDir "$workflow.projectDir", mode:'copy'

    output:
    tuple val("$workflow.projectDir/test_data/project_run.csv"), file('test_data/project_run.csv'), path('test_data')
    
    script:
    """
    wget -O test_data.tar.gz "https://www.dropbox.com/scl/fi/9dmjdg6f6kqpkgbxm87t7/test_data.tar.gz?rlkey=qdvk54r3dq7ib5t46pnuvjeka&st=b7fc5m5c&dl=0"
    tar -xzvf test_data.tar.gz
    """
}

process DOWNLOAD_REF_PANEL {
    publishDir "$workflow.outputDir", mode:'copy'

    output:
    path "1000G_hg38"
    
    script:
    """
    wget -O 1000G_hg38.zip "https://www.dropbox.com/scl/fi/d0xm2ht5scduspucmr2qb/1000G_hg38.zip?rlkey=niv01mpleiqjosvcljgguf6t9&st=9lxnh2tf&dl=0"
    unzip 1000G_hg38.zip
    """
}

process NUMBAT_PROC_RNA {
    publishDir "$workflow.outputDir", mode:'copy', pattern: '*_allele_counts.tsv.gz'

    input:
    tuple val(sample_id), path(rna_bam), path(rna_bai), path(rna_counts), path(rna_barcodes)
    each path(gmap)
    each path(snps)
    each path(panel)
    each path(p_script)

    output:
    tuple val(sample_id), path("count_alleles/$sample_id/${sample_id}_allele_counts.tsv.gz"), path(rna_counts), path("count_alleles/$sample_id/*.log")

    script:
    """
    MEM='$task.memory'
    mkdir -p count_alleles/$sample_id

    Rscript $p_script \
        --label $sample_id \
        --samples $sample_id \
        --bams $rna_bam \
        --barcodes $rna_barcodes \
        --outdir count_alleles/$sample_id \
        --gmap $gmap \
        --snpvcf $snps \
        --paneldir $panel \
        --ncores $task.cpus
    """
}

process NUMBAT_RUN_RNA {
    publishDir "$workflow.outputDir", mode:'copy'

    input:
    tuple val(sample_id), path(df_allele), path(rna_counts), path(logs)
    each path(script)

    output:
    tuple val(sample_id), path("numbat/$sample_id/*")

    script:
    """
    MEM='$task.memory'
    mkdir -p numbat/$sample_id

    Rscript $script \
        --countmat $rna_counts \
        --alleledf $df_allele \
        --ncores $task.cpus \
        --outdir numbat/$sample_id
    """
}

workflow {
    main:
    if (params.download_data) {
        DOWNLOAD_TEST_DATA()
        samples = DOWNLOAD_TEST_DATA.out.map { t -> t[1] }.splitCsv(header:true)
    } else {
        samples = channel.fromPath(params.sample_sheet)
                         .splitCsv(header:true)
    }

    gmap = channel.fromPath(params.eagle_gmap)
    if (params.numbat_matched_normal_phasing) {
        p_script = channel.fromPath('scripts/pileup_no_phase.R')
        matched_phasing(samples.map { row -> [row.sample, row.normal_bam, row.normal_bai] })
        numbat_rna(samples.map { row -> [row.sample, "${projectDir}/${row.rna_bam}", "${projectDir}/${row.rna_bai}", 
                                         "${projectDir}/${row.rna_counts}", "${projectDir}/${row.rna_barcodes}"] }, 
                    gmap, matched_phasing.out.phased_snps, p_script)
    } else {
        p_script = channel.fromPath('scripts/pileup_and_phase.R')
        snps = channel.fromPath(params.genome1k_snps)
        numbat_rna(samples.map { row -> [row.sample, "${projectDir}/${row.rna_bam}", "${projectDir}/${row.rna_bai}", 
                                         "${projectDir}/${row.rna_counts}", "${projectDir}/${row.rna_barcodes}"] }, 
                    gmap, snps, p_script)
    }

    // publish:
    // het_snps = matched_phasing.out.het_snps
    // phased_snps = matched_phasing.out.phased_snps
    // allele_counts = numbat_rna.out.allele_counts
    // numbat_rna_out = numbat_rna.out.numbat_rna_out
}

workflow matched_phasing {
    take:
    samples

    main:
    def ref_path = new File(params.ref)

    HAPLOTYPECALLER(
        samples, 
        ref_path.parent,
        ref_path.name,
        params.chroms
    )

    BCFTOOLS_HC(HAPLOTYPECALLER.out.groupTuple(size:22)) // TODO: fix hardcoded n chroms

    EAGLE2(
        BCFTOOLS_HC.out,
        params.eagle_ref_panel_dir,
        params.eagle_gmap,
        params.chroms
    )

    BCFTOOLS_EAGLE(EAGLE2.out.groupTuple(size:22))

    emit:
    het_snps = BCFTOOLS_HC.out
    phased_snps = BCFTOOLS_EAGLE.out
}

workflow numbat_rna {
    take: 
    samples
    gmap
    snps
    p_script

    main:
    script = channel.fromPath("scripts/run_numbat.R")
    if (params.numbat_matched_normal_phasing) {
        NUMBAT_PROC_RNA(samples, gmap, snps, null, p_script)
    } else {
        if (params.download_data) {
            panel = DOWNLOAD_REF_PANEL()
        } else {
            panel = channel.fromPath(params.eagle_ref_panel_dir)
        }
        NUMBAT_PROC_RNA(samples, gmap, snps, panel, p_script)
    }
    NUMBAT_RUN_RNA(NUMBAT_PROC_RNA.out, script)

    emit:
    allele_counts = NUMBAT_PROC_RNA.out
    numbat_rna_out = NUMBAT_RUN_RNA.out
}

// output {
//     het_snps {
//     }
//     phased_snps{
//     }
//     allele_counts{
//     }
//     numbat_rna_out{
//     }
// }
