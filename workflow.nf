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

workflow {
    main:
    samples = channel.fromPath(params.sample_sheet)
                     .splitCsv(header:true)

    if (params.numbat_matched_normal_phasing) {
        matched_phasing(samples.map { row -> [row.sample, row.normal_bam, row.normal_bai] })
    }

    publish: 
    het_snps = matched_phasing.out.het_snps
    phased_snps = matched_phasing.out.phased_snps
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

output {
    het_snps {
    }
    phased_snps{
    }
}