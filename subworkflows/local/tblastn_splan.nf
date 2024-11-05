//
// Buid BlastDB and align with related genomes
//

include { BLAST_MAKEBLASTDB } from '../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_TBLASTN } from '../../modules/nf-core/blast/tblastn/main'
include { EXONERATE }  from '../../modules/local/exonerate'
include { SEQTK_EXTRACT } from '../../modules/local/seqtk'
include { SPALN_MAKEINDEX } from '../../modules/local/spaln/makeindex'
include { SPALN_ALIGN } from '../../modules/local/spaln/align'
include { SPALN_MERGE } from '../../modules/local/spaln/merge'
include { HELPER_SPALNTOEVM } from '../../modules/local/helper/spalntoevm'
include { HELPER_SPALNTOGMOD } from '../../modules/local/helper/spalntogmod'
include { HELPER_SPALNTOTRAINING } from '../../modules/local/helper/spalntotraining'
include { FILTER_PROTEIN } from '../../modules/local/filter_protein'



workflow TBLASTN_EXONERATE{

    take:
    assembly
    ch_ref_genomes
    protein_identity


    main:

    ch_ref_genomes.splitCsv( header:true)
                  .map{ create_target_channel(it) }
		  .map{meta, assembly, proteins, gtf -> tuple(meta, proteins)  }
		  .set{ proteins_targets }

    BLAST_MAKEBLASTDB(assembly)

    SPALN_MAKEINDEX(assembly)  

    BLAST_TBLASTN(proteins_targets, BLAST_MAKEBLASTDB.out.db)

    tblastn_results = BLAST_TBLASTN.out.txt.join(proteins_targets)
    
    SEQTK_EXTRACT(tblastn_results)

    protein_chunks = SEQTK_EXTRACT.out.tblastn_proteins.splitFasta(by: params.nproteins, file: true, elem: [1])
       
    SPALN_ALIGN(
          SPALN_MAKEINDEX.out.spaln_index,
          protein_chunks,
          params.spaln_q,
          params.spaln_options)

     
    SPALN_MERGE(
          SPALN_MAKEINDEX.out.spaln_index,
          SPALN_ALIGN.out.align.groupTuple().map { m,files -> tuple(m,files.flatten()) },
          protein_identity)

    HELPER_SPALNTOGMOD(SPALN_MERGE.out.gff)
       
    HELPER_SPALNTOEVM(SPALN_MERGE.out.gff)
       
    HELPER_SPALNTOTRAINING(SPALN_MERGE.out.gff)
       
    emit:
    gff = SPALN_MERGE.out.gff
    gff_training = HELPER_SPALNTOTRAINING.out.gff
    evm = HELPER_SPALNTOEVM.out.gff

}




def create_target_channel(LinkedHashMap row) {

    //species, assembly, protiens, gtf
    def meta = [:]
    meta.id           = row.species
    
    return [ meta, file(row.assembly),file(row.protiens), file(row.gtf) ]
    
}

