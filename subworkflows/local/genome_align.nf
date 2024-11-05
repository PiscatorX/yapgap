//
// Align genomes and map annotations
//

include { FASTASPLITTER } from '../../modules/local/fastasplitter'
include { SATSUMA2_SATSUMASYNTENY2 } from '../../modules/local/satsuma2/satsumasynteny2'
include { KRAKEN } from '../../modules/local/kraken'
include { HELPER_KRAKEN2GFF as SATSUMA_KRAKEN2GFF } from '../../modules/local/helper/kraken2gff'
include { GAAS_FASTACLEANER } from '../../modules/local/gaas/fastacleaner'
include { HELPER_GTF2HINTS as SATSUMA_GTF2HINTS } from '../../modules/local/helper/gtf2hints'
include { GAAS_FASTAFILTERBYSIZE } from '../../modules/local/gaas/fastafilterbysize'

workflow GENOME_ALIGN{

    take:
    assembly 
    ch_ref_genomes
    
    main:
    //array = [ meta, file(assembly), file(proteins), file(gtf)
    ch_ref_genomes.splitCsv( header:true)
                  .map{ create_target_channel(it) }
    		  .set { targets }
		      
    assembly_targets = targets.map{meta, assembly, proteins, gtf -> tuple(meta, assembly) }
    
    GAAS_FASTACLEANER(assembly_targets)
    
    GAAS_FASTAFILTERBYSIZE(GAAS_FASTACLEANER.out.fasta, params.min_contig_size)

    // To DO: Replace with Nextflow SplitFasta
    FASTASPLITTER( assembly, params.npart_size )
   
    targets_clean = GAAS_FASTAFILTERBYSIZE.out.fasta.join(assembly_targets).map{ meta, filt_assembly, assembly -> tuple(meta, filt_assembly) }

    satsuma_input =  FASTASPLITTER.out.chunks.flatten()
                                             .map{it ->  [ ch_meta(it.getBaseName()), it ] }
			                     .combine(targets_clean)
			    
    SATSUMA2_SATSUMASYNTENY2(satsuma_input)

    grouped_chains = SATSUMA2_SATSUMASYNTENY2.out.satsuma_chain.groupTuple(by: [0,1,2,3])
    
    KRAKEN(grouped_chains)
    
    SATSUMA_KRAKEN2GFF(KRAKEN.out.gtf)
    
    SATSUMA_GTF2HINTS(KRAKEN.out.gtf,params.pri_trans)
    
    emit:
    versions = SATSUMA2_SATSUMASYNTENY2.out.versions
    gff = SATSUMA_KRAKEN2GFF.out.gff
    hints = SATSUMA_GTF2HINTS.out.gff


}


def ch_meta(basename){

    def meta = [:]
    meta.id  = basename
    
    return meta

}



def get_meta(f){

    def meta = [:]
    meta.id     = file(f).getSimpleName()
    
    return meta

}




def create_file_channel(f){

    def meta = [:]
    meta.id           = file(f).getSimpleName()

    def array = [ meta, file(f) ]

    return array

}



def create_target_channel(LinkedHashMap row) {

    //species, assembly, protiens, gtf
    def meta = [:]
    meta.id           = row.species
    
    return [ meta, file(row.assembly), file(row.protiens), file(row.gtf) ]
}




