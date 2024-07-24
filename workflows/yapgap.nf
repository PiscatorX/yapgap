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
log.info paramsSummaryLog(workflow)

// Create a new channel of metadata from a sample sheet
// NB: `input` corresponds to `params.input` and associated sample sheet schema
ch_input = Channel.fromSamplesheet("input")

WorkflowYAPGAP.initialise(params, log)




//def checkPathParamList = [ params.multiqc_config, params.assembly ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.assembly) { ch_assembly = file(params.assembly, checkIfExists: true) } else { exit 1, 'No assembly specified!' }

// Set relevant input channels
if (params.proteins) { ch_proteins = file(params.proteins, checkIfExists: true) } else { ch_proteins = Channel.empty() }
if (params.proteins_targeted) { ch_proteins_targeted = file(params.proteins_targeted, checkIfExists: true) } else { ch_proteins_targeted = Channel.empty() }
if (params.transcripts) { ch_t = file(params.transcripts) } else { ch_transcripts = Channel.empty() }
if (params.rnaseq_samples) { ch_rna_samplesheet = file(params.rnaseq_samples, checkIfExists: true) } else { ch_rna_samplesheet = Channel.empty() }
if (params.rm_lib) { ch_repeats = Channel.fromPath(file(params.rm_lib, checkIfExists: true)) } else { ch_repeats = Channel.empty()  }
if (params.aug_config_dir) { ch_aug_config_folder = file(params.aug_config_dir, checkIfExists: true) } else { ch_aug_config_folder = Channel.from(params.aug_config_container) }
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
ch_rfam_cm = file("${workflow.projectDir}/assets/rfam/14.2/Rfam.cm.gz", checkIfExists: true)
ch_rfam_family = file("${workflow.projectDir}/assets/rfam/14.2/family.txt.gz", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTASPLITTER } from '../modules/local/fastasplitter'
include { ASSEMBLY_PREPROCESS } from '../subworkflows/local/assembly_preprocess'
include { REPEATMASKER } from '../subworkflows/local/repeatmasker'
include { SPALN_ALIGN_PROTEIN ; SPALN_ALIGN_PROTEIN as SPALN_ALIGN_MODELS } from '../subworkflows/local/spaln_align_protein'
include { RNASEQ_ALIGN_HISAT } from '../subworkflows/local/rnaseq_align_hisat2'
include { MINIMAP_ALIGN_TRANSCRIPTS ; MINIMAP_ALIGN_TRANSCRIPTS as TRINITY_ALIGN_TRANSCRIPTS } from '../subworkflows/local/minimap_align_transcripts'
include { AUGUSTUS_PIPELINE } from '../subworkflows/local/augustus_pipeline'
include { PASA_PIPELINE } from '../subworkflows/local/pasa_pipeline'
include { GENOME_ALIGN } from '../subworkflows/local/genome_align'
include { EVM } from '../subworkflows/local/evm.nf'
include { FASTA_PREPROCESS as TRANSCRIPT_PREPROCESS } from '../subworkflows/local/fasta_preprocess'
include { BUSCO_QC } from '../subworkflows/local/busco_qc'
include { NCRNA } from '../subworkflows/local/ncrna'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// MODULE: Installed directly from nf-core
//
include { ASSEMBLYSCAN } from '../modules/nf-core/assemblyscan/main'
include { SAMTOOLS_MERGE } from '../modules/local/samtools/merge'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRINITY_GENOMEGUIDED } from '../modules/local/trinity/genomeguided'
include { AUGUSTUS_BAM2HINTS } from '../modules/local/augustus/bam2hints'
include { REPEATMODELER } from '../modules/local/repeatmodeler'
include { AUGUSTUS_STAGECONFIG } from '../modules/local/augustus/stageconfig'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


ch_assembly  = create_file_channel(params.assembly)


workflow YAPGAP {

     // Assembly preprocess and get some stats

     ASSEMBLY_PREPROCESS(ch_assembly)
     ASSEMBLYSCAN(ch_assembly)
   
     total_contig_length = ASSEMBLYSCAN.out.json
    	                  .map{ meta, json -> json }
    	                  .splitJson()
    	                  .first{ it.key == "total_contig_length" }
    	                  .map{ it.value }
			  
    fasta_npart_size =  total_contig_length.map{ it.intdiv(params.max_cpus) }
   
    FASTASPLITTER(ch_assembly, fasta_npart_size)

    fasta_chunks = FASTASPLITTER.out.meta.combine(FASTASPLITTER.out.chunks.flatten())
    /* ######################################################################### */

    Channel.fromSamplesheet("input").splitCsv( header: true )

// .map {
//             meta, fastq_1, fastq_2 ->
//                 if (!fastq_2) {
//                     return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
//                 } else {
//                     return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
//                 }
//         }
//         .groupTuple()
//         .map {
//             WorkflowRnaseq.validateInput(it)
//         }
//         .branch {
//             meta, fastqs ->
//                 single  : fastqs.size() == 1
//                     return [ meta, fastqs.flatten() ]
//                 multiple: fastqs.size() > 1
//                     return [ meta, fastqs.flatten() ]
//         }
//         .set { ch_fastq }




//  ch_rna_samplesheet

// if (params.rnaseq_samples) {
//      RNASEQ_ALIGN(
//         ASSEMBLY_PREPROCESS.out.fasta.collect(),
//         ch_samplesheet
//      )

      
       // // 
       // // MODULE: Merge all BAM files
       // //
       // RNASEQ_ALIGN.out.bam.map{ meta, bam ->
       //  new_meta = [:]
       //  new_meta.id = meta.ref
       //  tuple(new_meta,bam)
       // }.groupTuple(by:[0])
       // .set{bam_mapped}
       
}


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     COMPLETION EMAIL AND SUMMARY
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.summary(workflow, params, log)
// }

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */



def create_file_channel(f) {

    def meta = [:]
    meta.id           = file(f).getSimpleName()

    def array = [ meta, file(f) ]

    return array

}







    // ch_empty_gff = Channel.fromPath(params.dummy_gff)
    // ch_versions = Channel.empty()
    // ch_hints = Channel.empty()
    // ch_repeats_lib = Channel.empty()
    // ch_proteins_gff = Channel.from([])
    // ch_transcripts_gff = Channel.from([])
    // ch_genes_gff = Channel.empty()
    // ch_transcripts = Channel.empty()
    // ch_genome_rm = Channel.empty()
    // ch_proteins_fa = Channel.empty()
    // ch_busco_qc = Channel.empty()

    // // SUBWORKFLOW: Repeat modelling if no repeats are provided
    // //
    // if (!params.rm_lib && !params.rm_species) {
    //    REPEATMODELER(
    //       ASSEMBLY_PREPROCESS.out.fasta
    //    )
    //    ch_repeats = REPEATMODELER.out.fasta.map {m,fasta -> fasta}
    // }

    // //
    // // MODULE: Repeatmask the genome; if a repeat species is provided, use that - else the repeats in FASTA format
    // if (params.rm_species) {
    //    REPEATMASKER(
    //       ASSEMBLY_PREPROCESS.out.fasta,
    //       ch_repeats,
    //       params.rm_species,
    //       ch_rm_db
    //    )
    //    ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
    //    ch_genome_rm = REPEATMASKER.out.fasta
    // } else {
    //    REPEATMASKER(
    //       ASSEMBLY_PREPROCESS.out.fasta,
    //       ch_repeats,
    //       false,
    //       ch_rm_db
    //    )
    //    ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
    //    ch_genome_rm = REPEATMASKER.out.fasta
    // }

    // //
    // // SUBWORKFLOW: Turn transcript inputs to channel
    // //
    // if (params.transcripts) {
    //    TRANSCRIPT_PREPROCESS(
    //       ch_t
    //    )
    //    ch_transcripts = ch_transcripts.mix(TRANSCRIPT_PREPROCESS.out.fasta)
    // }

    // //
    // // MODULE: Stage Augustus config dir to be editable
    // //
    // AUGUSTUS_STAGECONFIG(ch_aug_config_folder)
    // ch_aug_config_folder = AUGUSTUS_STAGECONFIG.out.config_dir

    // //
    // // SUBWORKFLOW: Validate and pre-process the assembly
    // //
    // ASSEMBLY_PREPROCESS(
    //     ch_genome
    // )
    // ch_versions = ch_versions.mix(ASSEMBLY_PREPROCESS.out.versions)

    // // 
    // // SUBWORKFLOW: Search for ncRNAs
    // //
    // if (params.ncrna) {
    //    NCRNA(
    //       ASSEMBLY_PREPROCESS.out.fasta,
    //       ch_rfam_cm,
    //       ch_rfam_family
    //    )
    // }

    // //
    // // SUBWORKFLOW: Align genomes and map annotations
    // //
    // if (params.references) {
    //    GENOME_ALIGN(
    //       ASSEMBLY_PREPROCESS.out.fasta,
    //       ch_ref_genomes
    //    )
    //    ch_versions = ch_versions.mix(GENOME_ALIGN.out.versions)
    //    ch_hints = ch_hints.mix(GENOME_ALIGN.out.hints)
    //    ch_genes_gff = ch_genes_gff.mix(GENOME_ALIGN.out.gff)
    // }          

    // //  

    // //
    // // SUBWORKFLOW: Align proteins from related organisms with SPALN
    // if (params.proteins) {
    //    SPALN_ALIGN_PROTEIN(
    //       ASSEMBLY_PREPROCESS.out.fasta,
    //       ch_proteins,
    //       params.spaln_protein_id
    //    )
    //    ch_versions = ch_versions.mix(SPALN_ALIGN_PROTEIN.out.versions)
    //    ch_hints = ch_hints.mix(SPALN_ALIGN_PROTEIN.out.hints)
    //    ch_proteins_gff = ch_proteins_gff.mix(SPALN_ALIGN_PROTEIN.out.evm)
    // } 

    // // 
    // // SUBWORKFLOW: Align species-specific proteins 
    // if (params.proteins_targeted) {
    //    SPALN_ALIGN_MODELS(
    //       ASSEMBLY_PREPROCESS.out.fasta,
    //       ch_proteins_targeted,
    //       params.spaln_protein_id_targeted
    //    )
    //    ch_versions = ch_versions.mix(SPALN_ALIGN_MODELS.out.versions)
    //    ch_hints = ch_hints.mix(SPALN_ALIGN_MODELS.out.hints)
    //    ch_genes_gff = ch_genes_gff.mix(SPALN_ALIGN_MODELS.out.gff)
    // }

    // //
    // // SUBWORKFLOW: Align RNAseq reads
    // //
    // 
    //    //
    //    // MODULE: Merge BAM files
    //    //
    //    SAMTOOLS_MERGE(
    //       bam_mapped
    //    )
    //    AUGUSTUS_BAM2HINTS(
    //       SAMTOOLS_MERGE.out.bam,
    //       params.pri_rnaseq
    //    )
    //    ch_hints = ch_hints.mix(AUGUSTUS_BAM2HINTS.out.gff)
    //    ch_versions = ch_versions.mix(RNASEQ_ALIGN.out.versions.first(),AUGUSTUS_BAM2HINTS.out.versions,SAMTOOLS_MERGE.out.versions)

    //    //
    //    // SUBWORKFLOW: Assemble transcripts using Trinity and align to genome
    //    //
    //    if (params.trinity) {
    //       TRINITY_GENOMEGUIDED(
    //          SAMTOOLS_MERGE.out.bam,
    //          params.max_intron_size
    //       )
    //       ch_transcripts = ch_transcripts.mix(TRINITY_GENOMEGUIDED.out.fasta)
    //       ch_versions = ch_versions.mix(TRINITY_GENOMEGUIDED.out.versions)
    //    }
    // }

    // //
    // // SUBWORKFLOW: Align transcripts to the genome
    // //

    // if (params.transcripts || params.trinity) {
    //    MINIMAP_ALIGN_TRANSCRIPTS(
    //       ASSEMBLY_PREPROCESS.out.fasta.collect(),
    //       ch_transcripts
    //    )
    //    ch_versions = ch_versions.mix(MINIMAP_ALIGN_TRANSCRIPTS.out.versions)
    //    ch_transcripts_gff = ch_transcripts_gff.mix(MINIMAP_ALIGN_TRANSCRIPTS.out.gff)
    //    ch_hints = ch_hints.mix(MINIMAP_ALIGN_TRANSCRIPTS.out.hints)
    // }

    // //
    // // SUBWORKFLOW: Assemble transcripts into gene models
    // //
    // if (params.pasa) {
    //     PASA_PIPELINE(
    //        ASSEMBLY_PREPROCESS.out.fasta,
    //        ch_transcripts.map { m,t -> t }.collectFile(name: "transcripts.merged.fa").map { it ->
    //           def mmeta = [:]
    //           mmeta.id = "merged"
    //           tuple(mmeta,it)
    //        }
    //     )
    //     ch_versions = ch_versions.mix(PASA_PIPELINE.out.versions)
    //     ch_genes_gff = ch_genes_gff.mix(PASA_PIPELINE.out.gff)
    // }

    // //
    // // SUBWORKFLOW: Train augustus prediction model
    // //
    // if (params.aug_training) {

    //    if (params.proteins_targeted) {
    //       AUGUSTUS_TRAINING(
    //          SPALN_ALIGN_MODELS.out.gff
    //       )
    //       ch_aug_config_final = AUGUSTUS_TRAINING.out.aug_config_folder             
    //    } else if (params.pasa) {
    //       AUGUSTUS_TRAINING(
    //          PASA_PIPELINE.out.gff_training
    //       )
    //       ch_aug_config_final = AUGUSTUS_TRAINING.out.aug_config_folder
    //    }

    // } else {
    //    ch_aug_config_final = ch_aug_config_folder
    // }

    // //
    // // SUBWORKFLOW: Predict gene models using AUGUSTUS
    // //
    // all_hints = ch_hints.unique().collectFile(name: 'hints.gff')

    // AUGUSTUS_PIPELINE(
    //    REPEATMASKER.out.fasta,
    //    all_hints,
    //    ch_aug_config_folder,
    //    ch_aug_extrinsic_cfg,
    // )
    // ch_versions = ch_versions.mix(AUGUSTUS_PIPELINE.out.versions)
    // ch_genes_gff = ch_genes_gff.mix(AUGUSTUS_PIPELINE.out.gff)
    // ch_proteins_fa = ch_proteins_fa.mix(AUGUSTUS_PIPELINE.out.proteins)

    // //
    // // SUBWORKFLOW: Consensus gene building with EVM
    // //
    // if (params.evm) {
    //    EVM(
    //       ch_genome_rm,
    //       ch_genes_gff.map{m,g -> g}.collectFile(name: 'genes.gff3'),
    //       ch_proteins_gff.map{m,p -> p}.mix(ch_empty_gff).collectFile(name: 'proteins.gff3'),
    //       ch_transcripts_gff.map{m,t ->t}.mix(ch_empty_gff).collectFile(name: 'transcripts.gff3'),
    //       ch_evm_weights
    //    )
    //    ch_proteins_fa = ch_proteins_fa.mix(EVM.out.proteins)
    // }

    // //
    // // SUBWORKFLOW: Check proteome completeness with BUSCO
    // //
    // if (params.busco_lineage) {
    //    BUSCO_QC(
    //       ch_proteins_fa,
    //       params.busco_lineage,
    //       params.busco_db_path
    //    )
    //    ch_busco_qc = BUSCO_QC.out.busco_summary
    // } 

    // //
    // // MODULE: Collect all software versions
    // //

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    // //
    // // MODULE: MultiQC
    // //
    // workflow_summary    = WorkflowGenomeannotator.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(ch_busco_qc.collect().ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
    // ch_versions    = ch_versions.mix(MULTIQC.out.versions)

