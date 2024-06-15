#!/usr/bin/env nextflow

include { TRNASCANSE } from './modules/local/trnascanse'
include { ASSEMBLYSCAN } from './modules/nf-core/assemblyscan/main'
include { INFERNAL_PRESS } from './modules/local/infernal/press'
include { INFERNAL_SEARCH } from './modules/local/infernal/search'
include { FASTASPLITTER } from './modules/local/fastasplitter'
include { HELPER_RFAMTOGFF } from './modules/local/helper/rfamtogff'
include { GUNZIP as GUNZIP_RFAM_CM; GUNZIP as GUNZIP_RFAM_FAMILY } from './modules/nf-core/gunzip/main'
include { TOGFF } from './modules/local/output2gff'
//include { GT_GFF3VALIDATOR } from '../modules/nf-core/gt/gff3validator/main'


rfam_cm_gz = file("${workflow.projectDir}/assets/rfam/14.2/Rfam.cm.gz", checkIfExists: true)
rfam_family_gz = file("${workflow.projectDir}/assets/rfam/14.2/family.txt.gz", checkIfExists: true)

//assembly_path = Channel.fromPath(params.assembly)
assembly = create_file_channel(params.assembly)


workflow NCRNA {

   take:
   assembly
   rfam_cm_gz
   rfam_family_gz

   main:
   assembly_basename = assembly[1].baseName
   
   ASSEMBLYSCAN(assembly)
   
   total_contig_length = ASSEMBLYSCAN.out.json
	                 .map{ meta, json -> json }
	                 .splitJson()
	                 .first{ it.key == "total_contig_length" }
	                 .map{ it.value }
   fasta_npart_size =  total_contig_length.map{ it.intdiv(params.max_cpus) }
   
   FASTASPLITTER(   
      assembly,
      fasta_npart_size
   )

   fasta_chunks = FASTASPLITTER.out.fasta_meta
   | combine(FASTASPLITTER.out.chunks.flatten())

   TRNASCANSE(fasta_chunks)

   TRNASCANSE.out.trnascan
   | map{ meta, output -> output }
   | flatten()  
   | branch{
            fos: it.getExtension() == "fos" 
            fpass_out: it.getExtension() == "fpass_out"
	    other: true
	    }
   | set{ trnascan_results }

   //  trnascan_results.fos
   // | collectFile ( name: "${assembly_basename}_fos.fasta" , storeDir: "${params.outdir}/tRNAscan/" ){ it.text }
   
   trnascan_results.fpass_out
   | collectFile(keepHeader: true, skip: 2 ){ it.text }
   | first{ it }
   | set{ trnascan_combined }

   TOGFF(FASTASPLITTER.out.fasta_meta, trnascan_combined, "tRNAscan") 

   // GUNZIP_RFAM_CM(
   //    create_file_channel(rfam_cm_gz)
   // )
   // GUNZIP_RFAM_FAMILY(
   //    create_file_channel(rfam_family_gz)
   // )
   
   // INFERNAL_PRESS(		
   // 	GUNZIP_RFAM_CM.out.gunzip
   // )
   
   //db_size_mb = total_contig_length.map{ 2*(it/1000000)}   

  //https://docs.rfam.org/en/latest/genome-annotation.html#understanding-infernal-output
   // INFERNAL_SEARCH(
   //    FASTASPLITTER.out.fasta_meta,
   //    db_size_mb,	
   //    FASTASPLITTER.out.chunks,
   //    INFERNAL_PRESS.out.cm_file,
   //    INFERNAL_PRESS.out.cm_index
   // )
   
   //HELPER_RFAMTOGFF(INFERNAL_SEARCH.out.tbl, GUNZIP_RFAM_FAMILY.out.gunzip.map{m,f -> f})

   // HELPER_RFAMTOGFF.out.gff
   // | collectFile ( name: "${assembly_basename}_infernal.gff" , storeDir: "${params.outdir}/Infernal/" )       
}


workflow{

    NCRNA(assembly, rfam_cm_gz, rfam_family_gz)

}





def create_file_channel(f) {
    def meta = [:]
    meta.id           = file(f).getSimpleName()

    def array = [ meta, file(f) ]

    return array
}





// HELPER_RFAMTOGFF.out.gff
// | view{ ">> $it" }
   