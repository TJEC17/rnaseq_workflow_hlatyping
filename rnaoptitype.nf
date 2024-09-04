
process  QC_trim{

    publishDir('QC_trim', mode : 'copy')

    input:
    tuple val(sampid), path(fq)

    output:
    tuple val(sampid), path("${sampid}*_val_{1,2}.fq.gz")

    script:
    """
    trim_galore --paired --q 30 --fastqc $fq
    """
}

process Optitype{

    publishDir('Optitype_out', mode : 'copy')

    input:
    tuple val(sampid), val(reads)

    output:
    path"*"

    """

    python $baseDir/OptiTypePipeline.py -i ${reads[0]} ${reads[1]} --dna -v -o $sampid

    """
}


process  STARalign{

    publishDir 'STAR_out', mode : 'copy'

    input:
    tuple val(sampid), path(trimmed)

    output:
    tuple val(sampid), path("*Aligned.out.bam")

    script:
    """
    STAR --runThreadN 2 --readFilesCommand zcat --readFilesIn $trimmed --genomeDir $baseDir/reference/ --outSAMtype BAM Unsorted --outFileNamePrefix $sampid
    """
}

process  counts{

    publishDir 'FeatureCounts_out', mode : 'copy'

    input:
    tuple val(sampid), path(aligned)

    output:
    path("*.txt")

    script:
    """
    featureCounts -T 4 -a $baseDir/reference/Homo_sapiens.GRCh38.111.gtf -o ${sampid}_readcounts.txt $aligned -p
    """
}

workflow {
    fast="$baseDir/fastq_data/*_{1,2}.fastq.gz"
    fast_ch=Channel.fromFilePairs(fast)
    QC_trim(fast_ch)
    Optitype(QC_trim.out)
    STARalign(QC_trim.out)
    counts(STARalign.out)
}

