#!/usr/bin/env nextflow

/*
========================================================================================
    TES-FC PIPELINE
========================================================================================
    A Nextflow pipeline for detecting treatment failure in malaria using genomic data
    analysis and machine learning.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

params.outdir = params.outdir ?: "results"

if (!params.site) {
    log.error "Site parameter is required! Please specify --site <SITENAME>"
    System.exit(1)
}

def check_file_exists(file_path, name) {
    if (!file(file_path).exists()) {
        log.error "Required ${name} file not found: ${file_path}"
        System.exit(1)
    }
}

check_file_exists("${params.raw_data_dir}/metadata_tes_${params.site}.csv", "metadata")
check_file_exists("${params.raw_data_dir}/genomic_site_${params.site}.csv", "genomic site")
check_file_exists("${params.madhito_amps_file}", "madhito amplicons")

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { PRINT_HELP } from './modules/local/help'

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process DATA_CLEANING {
    tag "Cleaning data for ${params.site}"
    label 'process_medium'
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path metadata
    path genomic_site
    path madhito_amps
    path tes_dirs

    output:
    path "genomic_updated_${params.site}.csv", emit: genomic_updated
    path "metadata_updated_${params.site}.csv", emit: metadata_updated

    script:
    """
    export NXF_SITE="${params.site}"
    export NXF_MAF_FILTER="${params.maf_filter}"
    export NXF_MIN_ALLELE_READ_COUNT="${params.min_allele_read_count}"
    export NXF_CUM_CURVE_THRESHOLD="${params.cum_curve_threshold}"
    export NXF_OPTIM_AMPSET="${params.optim_ampset}"
    
    data_cleaning.R
    """
}

process COI_CALCULATION {
    tag "Calculating COI for ${params.site}"
    label 'process_low'
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path genomic_updated
    path metadata_updated

    output:
    path "genomic_updated_${params.site}.csv", emit: genomic_updated
    path "metadata_updated_${params.site}.csv", emit: metadata_updated

    script:
    """
    export NXF_SITE="${params.site}"
    coi_calculation.R
    """
}

process FNR_CALCULATION {
    tag "Calculating FNR for ${params.site}"
    label 'process_low'
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path genomic_updated
    path metadata_updated

    output:
    path "FNRs_${params.site}.RDS", emit: fnrs

    script:
    """
    export NXF_SITE="${params.site}"
    export NXF_FNR_02_03="${params.fnr_02_03}"
    export NXF_FNR_BELOW_02="${params.fnr_below_02}"
    fnr_calculation.R
    """
}

process CLONE_GENERATION {
    tag "Generating clones for ${params.site}"
    label 'process_medium'
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path genomic_updated
    path metadata_updated

    output:
    path "clones_genomic_data_${params.site}.csv", emit: clones_genomic

    script:
    """
    export NXF_SITE="${params.site}"
    export NXF_N_CLONES_TARGET="${params.n_clones_target}"
    clone_generation.R
    """
}

process MIX_GENERATION {
    tag "Generating mixtures for ${params.site}"
    label 'process_high'
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path clones_genomic
    path metadata_updated

    output:
    path "MIXES_METADATA_${params.site}.RDS", emit: mixes_metadata
    path "MIXES_GENOMIC_${params.site}.RDS", emit: mixes_genomic

    script:
    """
    export NXF_SITE="${params.site}"
    export NXF_INITIAL_SAMPLE_SIZE="${params.initial_sample_size}"
    export NXF_CLONE_SEED="${params.clone_seed}"
    mix_generation.R
    """
}

process PAIRS_GENERATION {
    tag "Generating pairs for ${params.site}"
    label 'process_medium'
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path metadata_updated
    path mixes_metadata
    path mixes_genomic

    output:
    path "PAIRS_METADATA_${params.site}.RDS", emit: pairs_metadata
    path "PAIRS_GENOMIC_${params.site}.RDS", emit: pairs_genomic

    script:
    """
    export NXF_SITE="${params.site}"
    pairs_generation.R
    """
}

process FEATURE_CALCULATION {
    tag "Calculating features for ${params.site}"
    label 'process_high'
    cpus { params.max_cpus }
    memory { params.max_memory }
    publishDir "${params.outdir}/results", mode: 'copy'

    input:
    path pairs_metadata
    path pairs_genomic
    path genomic_updated
    path metadata_updated

    output:
    path "${params.site}_TRAINING_DATA.csv", emit: training_data
    path "${params.site}_REAL_DATA.csv", emit: real_data

    script:
    """
    export NXF_SITE="${params.site}"
    export NXF_CHUNK_SIZE="${params.chunk_size}"
    export NXF_NUM_CORES="${task.cpus}"
    feature_calculation.R
    """
}

process TARGET_AUGMENTATION {
    tag "Augmenting training data for ${params.site}"
    label 'process_high'
    cpus { params.max_cpus }
    memory { params.max_memory }
    publishDir "${params.outdir}/results", mode: 'copy'

    input:
    path fnrs
    path pairs_metadata
    path pairs_genomic
    path clones_genomic
    path training_data

    output:
    path "${params.site}_TRAINING_DATA.csv", emit: training_data

    script:
    """
    export NXF_SITE="${params.site}"
    export NXF_CLONE_SEED="${params.clone_seed}"
    export NXF_PAIRS_SEED="${params.pairs_seed}"
    export NXF_CHUNK_SIZE="${params.chunk_size}"
    export NXF_NUM_CORES="${task.cpus}"
    target_augmentation.R
    """
}

process ML_MODELING {
    tag "Modeling outcomes for ${params.site}"
    label 'process_high'
    publishDir "${params.outdir}/results", mode: 'copy'

    input:
    path training_data
    path real_data

    output:
    path "GLMNET_model_${params.site}.RDS", emit: model
    path "${params.site}_REAL_DATA_PREDICTIONS.csv", emit: predictions
    path "training_Results_GLMNET_${params.site}.csv", emit: results
    path "${params.site}_model_results_GLMNET.png", emit: plots
    path "feat_importance_${params.site}.png", emit: importance

    script:
    """
    export NXF_SITE="${params.site}"
    export NXF_TRAIN_TEST_SPLIT="${params.train_test_split}"
    export NXF_CV_FOLDS="${params.cv_folds}"
    ml_modeling.R
    """
}

/*
========================================================================================
    NAMED WORKFLOW
========================================================================================
*/

workflow SMALARIA_FAILURE_DETECTION {

    main:

    if (params.help) {
        PRINT_HELP()
        return
    }

    log.info """
    ====================================
    TES-FC PIPELINE
    ====================================
    Site: ${params.site}
    Input directory: ${params.raw_data_dir}
    Output directory: ${params.outdir}
    ====================================
    """.stripIndent()

    metadata_ch      = Channel.fromPath("${params.raw_data_dir}/metadata_tes_${params.site}.csv")
    genomic_site_ch  = Channel.fromPath("${params.raw_data_dir}/genomic_site_${params.site}.csv")
    madhito_amps_ch  = Channel.fromPath("${params.madhito_amps_file}")
    tes_dirs_ch      = Channel.fromPath("${params.raw_data_dir}/*_RESULTS*", type: 'dir').filter { it.exists() }

    DATA_CLEANING(metadata_ch, genomic_site_ch, madhito_amps_ch, tes_dirs_ch.collect())
    COI_CALCULATION(DATA_CLEANING.out.genomic_updated, DATA_CLEANING.out.metadata_updated)
    FNR_CALCULATION(COI_CALCULATION.out.genomic_updated, COI_CALCULATION.out.metadata_updated)
    CLONE_GENERATION(COI_CALCULATION.out.genomic_updated, COI_CALCULATION.out.metadata_updated)
    MIX_GENERATION(CLONE_GENERATION.out.clones_genomic, COI_CALCULATION.out.metadata_updated)
    PAIRS_GENERATION(COI_CALCULATION.out.metadata_updated, MIX_GENERATION.out.mixes_metadata, MIX_GENERATION.out.mixes_genomic)
    FEATURE_CALCULATION(PAIRS_GENERATION.out.pairs_metadata, PAIRS_GENERATION.out.pairs_genomic, COI_CALCULATION.out.genomic_updated, COI_CALCULATION.out.metadata_updated)
    TARGET_AUGMENTATION(FNR_CALCULATION.out.fnrs, PAIRS_GENERATION.out.pairs_metadata, PAIRS_GENERATION.out.pairs_genomic, CLONE_GENERATION.out.clones_genomic, FEATURE_CALCULATION.out.training_data)
    ML_MODELING(TARGET_AUGMENTATION.out.training_data, FEATURE_CALCULATION.out.real_data)

    emit:
    model       = ML_MODELING.out.model
    predictions = ML_MODELING.out.predictions
    results     = ML_MODELING.out.results
}

/*
========================================================================================
    EXECUTE
========================================================================================
*/

workflow {
    SMALARIA_FAILURE_DETECTION()
}

/*
========================================================================================
    COMPLETION SUMMARY
========================================================================================
*/

workflow.onComplete {
    log.info """
    ====================================
    PIPELINE EXECUTION SUMMARY
    ====================================
    Site: ${params.site}
    Completed at: ${workflow.complete}
    Duration: ${workflow.duration}
    Success: ${workflow.success}
    Work dir: ${workflow.workDir}
    Results: ${params.outdir}
    ====================================
    """.stripIndent()
}

