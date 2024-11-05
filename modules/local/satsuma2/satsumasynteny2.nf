process SATSUMA2_SATSUMASYNTENY2 {

    tag "${meta_q.id} | ${meta_t.id}"
    label 'process_medium'    
    label 'satsuma'
    
    
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using this version of Satsuma2. Please use docker or singularity containers."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/satsuma2:20161123--h7d875b9_3':
        'quay.io/biocontainers/satsuma2:20161123--h7d875b9_3' }"

    input:
    tuple val(meta_q), path(query), val(meta_t), path(target)


    output:
    tuple val(meta_q), path("$meta_t.id"), emit: satsuma_folder
    path("$satsuma_chain"),           	 emit: satsuma_chain
    path "versions.yml",	         emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_t.id}"
    satsuma_chain = meta_q.id + "_" + meta_t.id + "_satsuma.chained_out"
    
    // Satsuma is extremely chatty, need to redirect logs to /dev/null
    
    """

    export SATSUMA2_PATH=/usr/local/bin

    SatsumaSynteny2 \
       ${args} -q $query \
       -t $target \
       -threads ${task.cpus} \
       -o ${prefix}  2>&1 >/dev/null

    cp ${prefix}/satsuma_summary.chained.out $satsuma_chain

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        satsuma2: 20161123
    END_VERSIONS
    """
}
