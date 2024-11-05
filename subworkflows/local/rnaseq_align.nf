include { FASTQ_ALIGN_HISAT2 } from '../../subworkflows/local/fastq_align_hisat2/main'
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge'


workflow RNASEQ_ALIGN{

    take:
    ch_filtered_reads
    fai
    star_index
    hisat2_index
    ch_assembly
  

    main:
    
    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    ch_genome_bam_merged          = Channel.empty() 
    ch_genome_bam_merged_csi      = Channel.empty()
    ch_versions 		  = Channel.empty()
    
    
   if (params.aligner == 'hisat2') {
      //
      // SUBWORKFLOW: Alignment with HISAT2
      //
      ch_hisat2_multiqc = Channel.empty()
      FASTQ_ALIGN_HISAT2(ch_filtered_reads,
			    hisat2_index,
			    ch_assembly)

     ch_genome_bam        = FASTQ_ALIGN_HISAT2.out.bam
     ch_genome_bam_index  = FASTQ_ALIGN_HISAT2.out.bai
     ch_samtools_stats    = FASTQ_ALIGN_HISAT2.out.stats
     ch_samtools_flagstat = FASTQ_ALIGN_HISAT2.out.flagstat
     ch_samtools_idxstats = FASTQ_ALIGN_HISAT2.out.idxstats
     ch_hisat2_multiqc    = FASTQ_ALIGN_HISAT2.out.summary

     if (params.bam_csi_index) {

	 ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.csi

     }
     
    
     ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

     }else{

     if (params.aligner == 'star') {

	ALIGN_STAR (
	    ch_strand_inferred_filtered_fastq,
	    PREPARE_GENOME.out.star_index.map { [ [:], it ] },
	    PREPARE_GENOME.out.gtf.map { [ [:], it ] },
	    params.star_ignore_sjdbgtf,
	    '',
	    params.seq_center ?: '',
	    is_aws_igenome,
	    PREPARE_GENOME.out.fasta.map { [ [:], it ] }
	)

	ch_genome_bam        = ALIGN_STAR.out.bam
	ch_genome_bam_index  = ALIGN_STAR.out.bai
	ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
	ch_samtools_stats    = ALIGN_STAR.out.stats
	ch_samtools_flagstat = ALIGN_STAR.out.flagstat
	ch_samtools_idxstats = ALIGN_STAR.out.idxstats
	ch_star_multiqc      = ALIGN_STAR.out.log_final

	if (params.bam_csi_index) {
	    ch_genome_bam_index = ALIGN_STAR.out.csi
	}

	ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)

     }

    }

     bam_files = ch_genome_bam.map{ it[1]}.collect().map{ collection_meta(it) }

     SAMTOOLS_MERGE(bam_files, ch_assembly, fai.collect().map{[[:], it ]} )
     
     ch_genome_bam_merged = SAMTOOLS_MERGE.out.bam       

     ch_genome_bam_merged_csi = SAMTOOLS_MERGE.out.csi

    
     emit:
     ch_star_multiqc
     ch_genome_bam       
     ch_genome_bam_index 
     ch_samtools_stats   
     ch_samtools_flagstat
     ch_samtools_idxstats
     ch_hisat2_multiqc
     ch_genome_bam_merged        
     ch_genome_bam_merged_csi
     ch_versions   


}




def collection_meta(collection){

    def meta = [:]
    meta.id = file(params.assembly).getSimpleName()
    def array = [ meta, collection ]

    return array
}
