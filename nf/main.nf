process GeneratePaths {
    input:
    path samples

    output:
    path "${params.bam_paths}"

    script:
    """
    python3 ${baseDir}/scripts/generate_paths.py \\
        --samples $samples \\
        --output "${params.bam_paths}"
    """

}

process UnlinkPaths {
    input:
    path bam_paths

    output:
    path "${params.real_bam_paths}"

    script:
    """
    python3 ${baseDir}/scripts/unlink_paths.py \\
        --bam_paths $bam_paths \\
        --output "${params.real_bam_paths}"
    """
}

process CheckPaths {
    input:
    path real_bam_paths

    output:
    path "${params.missing_paths}"

    script:
    """
    python3 ${baseDir}/scripts/check_paths.py \\
        --bam_paths $real_bam_paths \\
        --output "${params.missing_paths}"
    """
}

workflow {
    samples_ch = file(params.samples)

    GeneratePaths(samples_ch)
    UnlinkPaths(GeneratePaths.out)
    CheckPaths(UnlinkPaths.out)
}
