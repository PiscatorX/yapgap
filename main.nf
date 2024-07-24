#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

//params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

//WorkflowMain.initialise(workflow, params, log)


include { YAPGAP } from './workflows/yapgap'






workflow YAPGAP_WORKFLOW {

     YAPGAP()
}



//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {

    YAPGAP_WORKFLOW()
}

