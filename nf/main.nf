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

process GenerateClinicalPaths {
    input:
    path samples
    path id_mapping
    val dmp_dir
    val mirror_bam_dir
    val mirror_access_bam_dir

    path dmp_key_path, name: "dmp_key.txt"
    path access_key_path, name: "access_key.txt"

    output:
    path "${params.clinical_bams}"

    script:
    """
    python3 ${baseDir}/scripts/generate_clinical_paths.py \\
        --samples $samples \\
        --id_mapping $id_mapping \\
        --dmp_dir $dmp_dir \\
        --mirror_bam_dir $mirror_bam_dir \\
        --mirror_access_bam_dir $mirror_access_bam_dir \\
        --dmp_key_path dmp_key.txt \\
        --access_key_path access_key.txt \\
        --output "${params.clinical_bams}"
    """


}

workflow {
    samples_ch = file(params.samples)
    id_ch = file(params.id_mapping)

    dmp_dir = params.dmp_dir
    mirror_bam_dir = params.mirror_bam_dir
    mirror_access_bam_dir = params.mirror_access_bam_dir
    dmp_key_path = file(params.dmp_key_path)
    access_key_path = file(params.access_key_path)

    GeneratePaths(samples_ch)
    UnlinkPaths(GeneratePaths.out)
    CheckPaths(UnlinkPaths.out)

    GenerateClinicalPaths(
        samples = samples_ch, 
        id_mapping = id_ch, 
        dmp_dir = dmp_dir,
        mirror_bam_dir = mirror_bam_dir,
        mirror_access_bam_dir = mirror_access_bam_dir,
        dmp_key_path = dmp_key_path,
        access_key_path = access_key_path)

}
