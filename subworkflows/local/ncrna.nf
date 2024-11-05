#!/usr/bin/env nextflow

include { TRNASCANSE } from '../../modules/local/trnascanse'
include { ASSEMBLYSCAN } from '../../modules/nf-core/assemblyscan/main'
include { INFERNAL_PRESS } from '../../modules/local/infernal/press'
include { INFERNAL_SEARCH } from '../../modules/local/infernal/search'
include { HELPER_RFAMTOGFF } from '../../modules/local/helper/rfamtogff'
include { GUNZIP as GUNZIP_RFAM_CM; GUNZIP as GUNZIP_RFAM_FAMILY } from '../../modules/nf-core/gunzip/main'
include { TOGFF } from '../../modules/local/output2gff'
//include { GT_GFF3VALIDATOR } from '../../modules/nf-core/gt/gff3validator/main'



workflow NCRNA{

   take:
   ch_assembly
   total_contig_length
   fasta_chunks
   rfam_cm_gz
   rfam_family_gz
   
   main:

   GUNZIP_RFAM_CM(rfam_cm_gz)
   GUNZIP_RFAM_FAMILY(rfam_family_gz)


   INFERNAL_PRESS(GUNZIP_RFAM_CM.out.gunzip)
   
   db_size_mb = total_contig_length.map{ 2*(it/1000000)}   

  //https://docs.rfam.org/en/latest/genome-annotation.html#understanding-infernal-output
   INFERNAL_SEARCH(db_size_mb,
		   fasta_chunks,
    		   INFERNAL_PRESS.out.cm_index,
    		   INFERNAL_PRESS.out.cm_file)

   HELPER_RFAMTOGFF(INFERNAL_SEARCH.out.tbl, GUNZIP_RFAM_FAMILY.out.gunzip.map{m,f -> f})

   assembly_basename = file(params.assembly).baseName
   
   HELPER_RFAMTOGFF.out.gff | collectFile ( name: "${assembly_basename}_infernal.gff" , storeDir: "${params.outdir}/Infernal/" )

 
   // TRNASCAN
   
   TRNASCANSE(fasta_chunks)

   TRNASCANSE.out.trnascan
   | map{ meta, output -> output }
   | flatten()  
   branch{
            fos: it.getExtension() == "fos" 
            fpass_out: it.getExtension() == "fpass_out"
   	    other: true
   	    }
   | set{ trnascan_results }

   
   trnascan_results.fos | collectFile ( name: "${assembly_basename}_fos.fasta" , storeDir: "${params.outdir}/tRNAscan/" ){ it.text }
   
   trnascan_results.fpass_out | collectFile(name: "${assembly_basename}.fpass_out", storeDir: "${params.outdir}/tRNAscan/", keepHeader: true, skip: 2 ){ it.text }
		              | set{ trnascan_combined }
			      
   //not tested
   TOGFF(ch_assembly, trnascan_combined, "tRNAscan") 

 
}



