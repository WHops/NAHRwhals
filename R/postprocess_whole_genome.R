
#' @export
tsv_to_bed_regional_dominance <- function(nahrwhals_res_tsv, nahrwhals_bed){
    postprocess_main_bash = system.file('extdata', 'scripts', 'scripts_postprocess', 'postprocess_transform_all.sh', package='nahrwhals')

    outdir = dirname(nahrwhals_bed)
    SCRIPTDIR_BASE = dirname(postprocess_main_bash)
    # Run the bash script
    command = paste0('bash ', postprocess_main_bash, ' ', nahrwhals_res_tsv, ' ', nahrwhals_bed, ' ', outdir, ' ', SCRIPTDIR_BASE, ' >/dev/null 2>&1')

    run_silent(command)
}

