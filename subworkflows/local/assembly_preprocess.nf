
//
// Clean and filter assembly
//

include { GAAS_FASTACLEANER } from '../../modules/local/gaas/fastacleaner'
include { GAAS_FASTASTATISTICS } from '../../modules/local/gaas/fastastatistics'
include { GAAS_FASTAFILTERBYSIZE as GAAS_ASSEMBLYFILTERBYSIZE } from '../../modules/local/gaas/fastafilterbysize'

workflow ASSEMBLY_PREPROCESS {

    take:
    assembly
    
    main:   	
    GAAS_ASSEMBLYFILTERBYSIZE(assembly, params.min_contig_size)    
    GAAS_FASTACLEANER(GAAS_ASSEMBLYFILTERBYSIZE.out.fasta)
    GAAS_FASTASTATISTICS(GAAS_FASTACLEANER.out.fasta)

    emit:
    fasta = GAAS_FASTACLEANER.out.fasta
    stats = GAAS_FASTASTATISTICS.out.stats
    versions = GAAS_ASSEMBLYFILTERBYSIZE.out.versions
}

