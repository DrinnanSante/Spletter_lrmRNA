
nextflow.enable.dsl=2
NXF_CONDA_ENABLED=true


process NanoPlot {
    tag "${sample_id}"

    conda "/home/drinnan/miniconda3/envs/spletter"

    publishDir "/home/drinnan/publish/${sample_id}/", mode: 'copy', overwrite: false, pattern: '*.html'


    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path(reads)
        path("*.html")
        

    script:
        """
        NanoPlot --fastq ${reads}
        """
}




/// Qualimap module: Does the QC for the alignment, report goes to the publish directory. Output has a pdf and a txt version, currently only taking the report
/// output is in a directory, fix to get only file is to use mv to move file out of directory to pwd
/// TODO find out if there is a speed up, understand the limits a little better on real data

process qualimap {
    tag "${sample_id}"

    publishDir "/home/drinnan/publish/${sample_id}/", mode: 'copy', overwrite: false, pattern: 'report.pdf', saveAs: {filename -> "${sample_id}_qualimap.pdf"}


    input: 
        tuple val(sample_id), path(bam)

    output:
        path("report.pdf")

    script:
        """
        /home/drinnan/qualimap/qualimap_v2.3/qualimap bamqc -outfile report.pdf -bam ${bam}
        mv marked_stats/report.pdf .
        """
}

/// Unmapped module: first extract all unmapped reads with samtools using view filtering for the 4 flag
/// then use bedtools to turn this into a fastq file for kraken analysis. Dont use samtools for this, it doesn't work, i dont know why
/// am putting all the files to their own folder so I can run kraken in serial
process unmapped_read_extraction {
    tag "${sample_id}"

    publishDir "/home/drinnan/publish/tokraken/", mode: 'copy', overwrite: false, pattern: 'unmapped.fq', saveAs: {filename -> "${sample_id}_kraken2_unmapped.fq"}

    input:
        tuple val(sample_id), path(bam)

    output:
        path("unmapped.fq")

    script:
    """
    samtools view -b -f 4 ${bam} > unmapped.bam
    bedtools bamtofastq -i unmapped.bam -fq unmapped.fq
    """
}

/// Kraken2 modlue: The only time I ran this it crashed my computer, cant use it until I find a way to use nextflow better. It ran all samples at once and the machine ran out of ram

process kraken2 {
    tag "${sample_id}"

    publishDir "/home/drinnan/publish/${sample_id}/", mode: 'copy', overwrite: false, pattern: 'contamreport.txt', saveAs: {filename -> "${sample_id}_kraken2.txt"}

    input:
        tuple val(sample_id), path(read)

    output:
        path("contamreport.txt")

    script:
    """
    /media/drinnan/genomes/Kraken_database/Kraken2/kraken2/kraken2 --db /media/drinnan/genomes/Kraken_database/Kraken2/database/ ${read}  --threads 30 --output kraken.txt --report contamreport.txt 
    """
}

/// This just takes way to long to use in practice

process kraken_uniq {
    tag "${sample_id}"

    publishDir "/home/drinnan/publish/${sample_id}/", mode: 'copy', overwrite: false, pattern: 'uniqreport.txt', saveAs: {filename -> "${sample_id}_kraken_uniq_unmapped.txt"}

    input:
        tuple val(sample_id), path(unmapped_fastq)

    output:
        path("uniqreport.txt")

    script:
    """
    /home/drinnan/Kraken/krakenuniq --db /media/drinnan/genomes/Kraken_database/ --threads 24  --output off --report-file uniqreport.txt bedtools.fq
    """
}
/*
/// Braken module: This currently doesn't work as you can't build the braken database. Author says they are working on it, check back later
process braken {
    tag "${sample_id}"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path()

    script:
    """
    samtools
    """
}
*/



/// BamIndex module: needed for other tools, don't want to carry it around unneeded modules

process bam_index {
    tag "${sample_id}"

    publishDir "/home/drinnan/publish/${sample_id}/", mode: 'copy', overwrite: false, pattern: 'marked.bam', saveAs: {filename -> "${sample_id}.bam"}
    publishDir "/home/drinnan/publish/${sample_id}/", mode: 'copy', overwrite: false, pattern: 'marked.bam.bai', saveAs: {filename -> "${sample_id}.bam.bai"}

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path(bam), path("marked.bam.bai")

    script:
        """
        samtools index ${bam}
        """
}


///TODO config file, sitag "${sample_id}"ngularity container, annotation, contamination check before mapping, contamination check on unmapped reads, somatic calling, SV calling, CNV calling?, repeat calling?, logic doc, somatic calling to find subcultures
///TODO logic doc, graph mappability, layout contamination checks, layout depthchecks, layout setup, layout F1 scores to show what can be found, show logical consitanty in controls 
workflow{
    reads = Channel.fromPath('/home/drinnan/projects/spletter_lrRNA/test_data/downsampled_10000.fastq')
              .map { file -> 
                  def sample_id = file.baseName.replaceAll(/_R\d+$/, '') 
                  [sample_id, file]
              }

    

    
    



    NanoPlot(reads)
   /*
    fastp(fastqc.out[0])
    
    
    bwa(fastp.out[0])
    mark_duplicates_picard(bwa.out)
    unmapped_read_extraction(mark_duplicates_picard.out[0])
    qualimap(mark_duplicates_picard.out[0])
    bam_index(mark_duplicates_picard.out[0])
    
    gridss(bam_index.out)
    gridss_cleanup(gridss.out)
    
    octopus_variant_call(bam_index.out)
    snpEff(octopus_variant_call.out)
    snpSift(snpEff.out)
    
    snp_call_check(octopus_variant_call.out) 
    */
    
}
    
    
