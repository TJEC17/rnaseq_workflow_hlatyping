
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


process seq2 {

    publishDir("seq2HLA_output", mode : 'copy')

    input:

    tuple val(sample_id), path(seq2)

    output:

    path "*"

    script:
    """ 

    python2 /usr/local/modules/easybuild/software/seq2HLA/2.3/bin/seq2HLA -1 ${seq2[0]} -2 ${seq2[1]} -r $sample_id

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
    seq2(QC_trim.out)
    STARalign(QC_trim.out)
    counts(STARalign.out)
}

