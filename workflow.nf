process HAPLOTYPECALLER {
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(ref_dir)
    val(ref_name)
    each chrom

    output:
    tuple val(sample_id), val(chrom)

    script:
    """
    REF=$ref_dir/$ref_name
    MEM='$task.memory'
    mkdir tmp
    BAM=\$(basename -s .bai $bai)
    mv $bam \$BAM

    gatk --java-options "-Xmx\${MEM%% *}g -XX:-UsePerfData" \
    HaplotypeCaller \
        --input \$BAM \
        --output ${sample_id}_${chrom}.vcf.gz \
        --reference \$REF \
        --native-pair-hmm-threads $task.cpus \
        --intervals $chrom \
        --tmp-dir tmp
    """
}

process CONCAT_SNPS {
    input:
    tuple val(sample_id), val(vcfs)

    output:
    tuple val(sample_id)

    script:
    """
    
    """
}

workflow {
    samples = channel.fromPath(params.sample_sheet)
                     .splitCsv(header:true)

    if (params.numbat_matched_normal_phasing) {
        matched_phasing(samples.map { row -> [row.sample, row.normal_bam, row.normal_bai] })
    }
}

workflow matched_phasing {
    take:
    samples

    main:
    def ref_path = new File(params.ref)
    HAPLOTYPECALLER(samples, 
                    ref_path.parent,
                    ref_path.name,
                    ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22'])
    CONCAT_SNPS(HAPLOTYPECALLER.out.groupTuple(size:22))

    emit:
    CONCAT_SNPS.out
}