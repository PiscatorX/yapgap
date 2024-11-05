/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
//log.info paramsSummaryLog(workflow)

WorkflowYAPGAP.initialise(params, log)

//def checkPathParamList = [ params.multiqc_config, params.assembly ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.assembly) { ch_assembly = file(params.assembly, checkIfExists: true) } else { exit 1, 'No assembly specified!' 
}
 

// Set relevant input channels

ch_rna_samplesheet = Channel.fromSamplesheet("rnaseq_samples")

if (params.proteins) { ch_proteins = file(params.proteins, checkIfExists: true) } else { ch_proteins = Channel.empty() }
if (params.busco_lineage) { ch_busco_lineage =  channel.fromPath(params.busco_lineage, checkIfExists: true) } else { ch_busco_lineage = Channel.empty() }
if (params.transcripts) { ch_t = file(params.transcripts) } else { ch_transcripts = Channel.empty() }
if (params.rm_lib) { ch_repeats = Channel.fromPath(file(params.rm_lib, checkIfExists: true)) } else { ch_repeats = Channel.empty()  }
if (params.references) { ch_ref_genomes = Channel.fromPath(params.references, checkIfExists: true)  } else { ch_ref_genomes = Channel.empty() }
if (params.rm_db)  { ch_rm_db = file(params.rm_db) } else { ch_rm_db = Channel.empty() }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_aug_extrinsic_cfg = params.aug_extrinsic_cfg ? Channel.from( file(params.aug_extrinsic_cfg, checkIfExists: true) ) : Channel.from( file("${workflow.projectDir}/assets/augustus/augustus_default.cfg"))
ch_evm_weights = Channel.from(file(params.evm_weights, checkIfExists: true))

ch_rfam_cm = Channel.fromPath(params.rfam_cm_gz, checkIfExists: true)
rfam_family_gz = Channel.fromPath(params.rfam_family_gz, checkIfExists: true)
ch_augustus_config_path = channel.fromPath(params.augustus_config_path, checkIfExists: true ) 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTASPLITTER } from '../modules/local/fastasplitter'
include { FASTP } from '../modules/nf-core/fastp/main'
include { ASSEMBLY_PREPROCESS } from '../subworkflows/local/assembly_preprocess'
include { REPEATMASKER } from '../subworkflows/local/repeatmasker'
include { SPALN_ALIGN_PROTEIN ; SPALN_ALIGN_PROTEIN as SPALN_ALIGN_MODELS } from '../subworkflows/local/spaln_align_protein'
include { MINIMAP_ALIGN_TRANSCRIPTS ; MINIMAP_ALIGN_TRANSCRIPTS as TRINITY_ALIGN_TRANSCRIPTS } from '../subworkflows/local/minimap_align_transcripts'
include { PASA_PIPELINE } from '../subworkflows/local/pasa_pipeline'
include { GENOME_ALIGN } from '../subworkflows/local/genome_align'
include { EVM } from '../subworkflows/local/evm.nf'
include { FASTA_PREPROCESS as TRANSCRIPT_PREPROCESS } from '../subworkflows/local/fasta_preprocess'
include { BUSCO_QC } from '../subworkflows/local/busco_qc'
include { NCRNA } from '../subworkflows/local/ncrna'
include { RNASEQ_INIT } from '../subworkflows/local/rnaseq_init'
include { RNASEQ_ALIGN } from '../subworkflows/local/rnaseq_align'
include { HISAT2_ALIGN } from '../modules/local/hisat2/align/main'
include { TRANSDECODER } from  '../subworkflows/local/transdecoder'
include { BRAKER } from  '../modules/local/braker'
include { TBLASTN_SPALN } from  '../subworkflows/local/tblastn_spaln'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core
//

include { AGAT_CONVERTSPGXF2GXF } from '../modules/nf-core/agat/convertspgxf2gxf/main'
include { BUSCO_BUSCO as BUSCO_GENOME; BUSCO_BUSCO as BUSCO_PROTEIN } from '../modules/nf-core/busco/busco/main'
include { ASSEMBLYSCAN } from '../modules/nf-core/assemblyscan/main'
include { SAMTOOLS_MERGE } from '../modules/local/samtools/merge'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRINITY_GENOMEGUIDED } from '../modules/local/trinity/genomeguided'
include { REPEATMODELER } from '../modules/local/repeatmodeler'
include { CAT_FASTQ } from '../modules/nf-core/cat/fastq'
include { FASTQ_FASTQC_UMITOOLS_FASTP } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { FASTQ_ALIGN_HISAT2 } from '../subworkflows/local/fastq_align_hisat2/main'
include { STRINGTIE_STRINGTIE } from '../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE } from '../modules/nf-core/stringtie/merge/main'
include { INTERPROSCAN } from '../modules/nf-core/interproscan/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_assembly  = create_file_channel(params.assembly)

workflow YAPGAP {

    ch_proteins_fa = Channel.empty()
    ch_busco_qc = Channel.empty()
    ch_versions = Channel.empty()
    ch_hints = Channel.empty()
    ch_genes_gff = Channel.empty()
    ch_proteins_gff = Channel.empty()
    ch_transcripts_gff = Channel.empty()		    

    /* ############################### PREPROCESS ################################### */
    
    ASSEMBLY_PREPROCESS(ch_assembly)
    ch_versions = ch_versions.mix(ASSEMBLY_PREPROCESS.out.versions)
    ASSEMBLYSCAN(ch_assembly)

    /* ################################## NCRNA ####################################### */

    //
    // SUBWORKFLOW: Check proteome completeness with BUSCO
    //

    if (params.genome_busco){

        lineage = channel.fromPath(params.busco_lineage).map{ it.baseName }

        BUSCO_GENOME(ch_assembly, "genome", lineage, ch_busco_lineage, [])
	
    }
   
    if (params.ncrna){
     
    	total_contig_length = ASSEMBLYSCAN.out.json
    			    .map{ meta, json -> json }
    			    .splitJson()
    			    .first{ it.key == "total_contig_length" }
    			    .map{ it.value }

    	fasta_npart_size =  total_contig_length.map{ it.intdiv( params.npart_size ) }

    	FASTASPLITTER(ch_assembly, fasta_npart_size)

    	fasta_chunks = FASTASPLITTER.out.meta.combine(FASTASPLITTER.out.chunks.flatten())
	
    	NCRNA(ch_assembly, total_contig_length, fasta_chunks, create_file_channel(params.rfam_cm_gz), create_file_channel(params.rfam_family_gz))
	
    }
    

    /* ############################# GENOME ALIGMENTS ##################################### */

    
    if (params.genome_align){

        GENOME_ALIGN(ch_assembly, ch_ref_genomes)  

    }


   /* ############################ PROTEIN ALIGNMENTS ##################################### */
   
  
   if (params.proteins_targeted){

        TBLASTN_SPALN(ch_assembly, ch_ref_genomes, params.spaln_protein_id_targeted)

	ch_genes_gff = ch_genes_gff.mix(TBLASTN_SPALN.out.gff)
		
     }
    
   
    if (params.other_proteins) {
     
        SPALN_ALIGN_PROTEIN(ASSEMBLY_PREPROCESS.out.fasta,
	                    ch_proteins,
	                    params.spaln_protein_id_targeted)
       
       ch_versions = ch_versions.mix(SPALN_ALIGN_PROTEIN.out.versions)
       ch_hints = ch_hints.mix(SPALN_ALIGN_PROTEIN.out.hints)
       ch_proteins_gff = ch_proteins_gff.mix(SPALN_ALIGN_PROTEIN.out.evm)
     
    } 


   /* ############################# RNASEQ DATA ##################################### */

    
    //RNA ALIGNMENTS
    // https://github.com/nf-core/rnaseq

    if (params.rnaseq) {
    
        // Create a new channel of metadata from a sample sheet

        def prepareToolIndices  = []
    	prepareToolIndices << params.aligner

        RNASEQ_INIT(ch_assembly,
    		     params.star_index,
    		     params.salmon_index,
    		     params.hisat2_index,
    		     prepareToolIndices)

        fai = RNASEQ_INIT.out.fai
    	chrom_sizes = RNASEQ_INIT.out.chrom_sizes 
    	star_index = RNASEQ_INIT.out.star_index             // channel: path(star/index/)
    	hisat2_index = RNASEQ_INIT.out.hisat2_index           // channel: path(hisat2/index/)

    	ch_versions = ch_versions.mix(RNASEQ_INIT.out.versions)

        ch_rna_samplesheet.map {meta, fastq_1, fastq_2 ->
	
    		    if (!fastq_2) {
    				return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
    			    } else {
    				return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
    			    }
    		    }
    		    .groupTuple()
    		    .map {
    			WorkflowRnaseq.validateInput(it)
    		    }
    		    .branch {
    			meta, fastqs ->
    			    single  : fastqs.size() == 1
    				return [ meta, fastqs.flatten() ]
    			    multiple: fastqs.size() > 1
    				return [ meta, fastqs.flatten() ]
    		    }
    		   .set { ch_fastq }


         CAT_FASTQ(ch_fastq.multiple)
     	           .reads
     	           .mix(ch_fastq.single)
     	           .set{ ch_cat_fastq }


        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    
       // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    
        ch_filtered_reads      = Channel.empty()
        ch_fastqc_raw_multiqc  = Channel.empty()
        ch_fastqc_trim_multiqc = Channel.empty()
        ch_trim_log_multiqc    = Channel.empty()
        ch_trim_read_count     = Channel.empty()

        if (params.trimmer == 'trimgalore') {

    	    FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
    	        ch_cat_fastq,
		params.skip_fastqc || params.skip_qc,
		params.with_umi,
		params.skip_umi_extract,
		params.skip_trimming,
		params.umi_discard_read,
		params.min_trimmed_reads)

        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
	ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
	ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
	ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
	ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
	ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
	
    }
    else{
    
    if (params.trimmer == 'fastp') {
     	      FASTQ_FASTQC_UMITOOLS_FASTP (
     		  ch_cat_fastq,
     		  params.skip_fastqc || params.skip_qc,
     		  params.with_umi,
     		  params.skip_umi_extract,
     		  params.umi_discard_read,
     		  params.skip_trimming,
     		  [],
     		  params.save_trimmed,
     		  params.save_trimmed,
     		  params.min_trimmed_reads
     	      )
	      
     	      ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
     	      ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
     	      ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
     	      ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
     	      ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
     	      ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    	      }
	      
     }

     RNASEQ_ALIGN(ch_filtered_reads,
                  fai,            
    		  star_index,
    		  hisat2_index,
    		  ch_assembly)

     ch_genome_bam        = RNASEQ_ALIGN.out.ch_genome_bam       
     ch_genome_bam_index  = RNASEQ_ALIGN.out.ch_genome_bam_index
     ch_samtools_stats    = RNASEQ_ALIGN.out.ch_samtools_stats
     ch_samtools_flagstat = RNASEQ_ALIGN.out.ch_samtools_flagstat
     ch_samtools_idxstats = RNASEQ_ALIGN.out.ch_samtools_idxstats
     ch_hisat2_multiqc    = RNASEQ_ALIGN.out.ch_hisat2_multiqc
     ch_genome_bam_merged = RNASEQ_ALIGN.out.ch_genome_bam_merged       
     ch_genome_bam_merged_csi = RNASEQ_ALIGN.out.ch_genome_bam_merged_csi
     //ch_versions = ch_versions.mix(RNASEQ_ALIGN.out.versions)


     if (params.stringtie) {

        STRINGTIE_STRINGTIE(ch_genome_bam, [])
	
	ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

        STRINGTIE_MERGE(STRINGTIE_STRINGTIE.out.transcript_gtf.collect { it[1] }, [])
        
	ch_transcripts_gff = AGAT_CONVERTSPGXF2GXF(STRINGTIE_MERGE.out.gtf.map{ gtf -> [insert_meta(params.assembly), gtf ] })

        TRANSDECODER(STRINGTIE_MERGE.out.gtf.map{ [ insert_meta(params.assembly), it ] },
     		    ch_assembly,
       		    params.Pfam_hmm,
       		    create_file_channel(params.Uniprot_ref))
		  
        cds  = TRANSDECODER.out.cds
        pep  = TRANSDECODER.out.pep
        bed  = TRANSDECODER.out.bed
        ch_proteins_gff = ch_proteins_gff.mix(TRANSDECODER.out.gff3)
        ch_versions = ch_versions.mix(TRANSDECODER.out.versions)
            
   if (params.pasa) {
    
        PASA_PIPELINE(ASSEMBLY_PREPROCESS.out.fasta, TRANSDECODER.out.transcripts)
        ch_versions = ch_versions.mix(PASA_PIPELINE.out.versions)
        ch_genes_gff = ch_genes_gff.mix(PASA_PIPELINE.out.gff)

     }
     
   }
   
   if (params.braker){
     			    
       // BRAKER(ch_assembly,
       //        params.protein_db,
       // 	      ch_genome_bam_merged,
       // 	      ch_augustus_config_path,
       // 	      [])
       // braker_gff = Channel.fromPath("/home/andhlovu/DB_REF/braker/braker.gff3")
       //ch_genes_gff = ch_genes_gff.mix(braker_gff)

    }
   }
   
    // if (params.evm) {
    //    EVM(ch_assembly,
    //        ch_genes_gff.map{m,g -> g}.collectFile(name: 'genes.gff3'),
    //        ch_proteins_gff.map{m,p -> p}.mix(ch_empty_gff).collectFile(name: 'proteins.gff3'),
    //        ch_transcripts_gff.map{m,t ->t}.mix(ch_empty_gff).collectFile(name: 'transcripts.gff3'),
    //        ch_evm_weights)
	   
    //    ch_proteins_fa = ch_proteins_fa.mix(EVM.out.proteins)
    //}

}





def get_meta(id, ch){

    return id

}



def insert_meta(f){

    def meta = [:]
    meta.id  = file(f).getSimpleName()

    return meta

}



def add_meta_abs(ch){

    def meta = [:]
    meta.id = ch.toAbsolutePath().toString()
    def array = [ meta, ch ]

    return array

}



def create_file_channel(f){

    def meta = [:]
    meta.id           = file(f).getSimpleName()

    def array = [ meta, file(f) ]

    return array

}



def empty_file_channel(){

    def meta = [:]
    meta.id           = ''

    def array = [ meta, ]

    return array

}