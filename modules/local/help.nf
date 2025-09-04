process PRINT_HELP {
    tag 'help'
    
    script:
    """
    cat << 'EOF'
    
    ====================================
    SMALARIA FAILURE DETECTION PIPELINE
    ====================================
    
    A Nextflow pipeline for detecting treatment failure in malaria using genomic data
    analysis and machine learning.
    
    Usage:
        nextflow run main.nf --site <SITENAME> [options]
    
    Required parameters:
        --site <string>                 Site name (e.g., 'Zambezia')
    
    Input/output options:
        --raw_data_dir <path>           Path to raw data directory (default: 'data/raw')
        --outdir <path>                 Output directory (default: 'results')
        --madhito_amps_file <path>      Path to madhito amplicons file 
                                       (default: 'resources/madhito_20amps.csv')
    
    Data processing parameters:
        --maf_filter <float>            Minor allele frequency threshold (default: 0.01)
        --min_allele_read_count <int>   Minimum reads per allele (default: 10)
        --cum_curve_threshold <float>   Cumulative curve threshold (default: 0.999)
        --optim_ampset <boolean>        Use optimized amplicon set (default: false)
    
    Clone generation parameters:
        --initial_sample_size <int>     Initial sample size for mixtures (default: 200)
        --n_clones_target <int>         Target number of clones (default: 1000)
    
    Random seeds:
        --clone_seed <int>              Seed for clone sampling (default: 420)
        --pairs_seed <int>              Seed for pairs sampling (default: 42069)
    
    ML model parameters:
        --train_test_split <float>      Training/test split ratio (default: 0.7)
        --cv_folds <int>                Cross-validation folds (default: 10)
    
    FNR parameters:
        --fnr_02_03 <float>            FNR for strains 0.02-0.03 abundance (default: 0.11)
        --fnr_below_02 <float>         FNR for strains below 0.02 abundance (default: 0.16)
    
    Processing parameters:
        --chunk_size <int>              Processing chunk size (default: 1000)
        --max_memory <string>           Maximum memory per process (default: '128.GB')
        --max_cpus <int>                Maximum CPUs per process (default: 16)
        --max_time <string>             Maximum time per process (default: '240.h')
    
    Profiles:
        -profile <string>               Configuration profile to use
                                       Available: conda, mamba, docker, singularity, 
                                                 podman, shifter, charliecloud, test
    
    Examples:
        # Basic run with required parameters
        nextflow run main.nf --site Zambezia
        
        # Run with custom data directory and output
        nextflow run main.nf --site Zambezia --raw_data_dir /path/to/data --outdir /path/to/results
        
        # Run with optimized amplicon set and more clones
        nextflow run main.nf --site Zambezia --optim_ampset true --n_clones_target 2000
        
        # Run with Docker
        nextflow run main.nf --site Zambezia -profile docker
    
    Input data requirements:
        The pipeline expects the following files in the raw data directory:
        
        1. metadata_tes_<SITE>.csv       - Sample metadata with PairsID and time_point
        2. genomic_site_<SITE>.csv       - Site genomic reference data
        3. <DIR>_RESULTS*/               - TES result directories containing
           allele_data_global_max_0_filtered.csv files
        4. resources/madhito_20amps.csv  - Amplicon reference file
    
    Output files:
        Key results will be saved to the output directory:
        
        - <SITE>_REAL_DATA_PREDICTIONS.csv     - Final predictions for TES pairs
        - GLMNET_model_<SITE>.RDS              - Trained ML model
        - training_Results_GLMNET_<SITE>.csv   - Model performance metrics
        - <SITE>_model_results_GLMNET.png      - Performance plots
        - feat_importance_<SITE>.png           - Feature importance plot
    
    For more information:
        https://github.com/smalaria/smalaria-failure-detection
    
    ====================================
    
EOF
    """
}
