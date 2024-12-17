process GenerateBamPaths {
    input:
    path samples
    path id_mapping

    val bam_path_normal
    val bam_path_plasma_duplex
    val bam_path_plasma_simplex
    val bam_path_index_normal
    val bam_path_index_duplex
    val bam_path_index_simplex
    val maf_path
    val cna_path
    val sv_path

    val dmp_dir
    val mirror_bam_dir
    val mirror_access_bam_dir

    path dmp_key_path, name: "dmp_key.txt"
    path access_key_path, name: "access_key.txt"

    script:
    """
    python3 ${baseDir}/scripts/per_patient.py \\
        --samples $samples \\
        --id_mapping $id_mapping \\
        --bam_path_normal $bam_path_normal \\
        --bam_path_plasma_duplex $bam_path_plasma_duplex \\
        --bam_path_plasma_simplex $bam_path_plasma_simplex \\
        --bam_path_index_normal $bam_path_index_normal \\
        --bam_path_index_duplex $bam_path_index_duplex \\
        --bam_path_index_simplex $bam_path_index_simplex \\
        --maf_path $maf_path \\
        --cna_path $cna_path \\
        --sv_path $sv_path \\
        --dmp_dir $dmp_dir \\
        --mirror_bam_dir $mirror_bam_dir \\
        --mirror_access_bam_dir $mirror_access_bam_dir \\
        --dmp_key_path dmp_key.txt \\
        --access_key_path access_key.txt \\
    """

}


workflow {
    samples_ch = file(params.samples)
    id_ch = file(params.id_mapping)

    bam_path_normal = params.bam_path_normal
    bam_path_plasma_duplex = params.bam_path_plasma_duplex
    bam_path_plasma_simplex = params.bam_path_plasma_simplex
    bam_path_index_normal = params.bam_path_index_normal
    bam_path_index_duplex = params.bam_path_index_duplex
    bam_path_index_simplex = params.bam_path_index_simplex

    maf_path = params.maf_path
    cna_path = params.cna_path
    sv_path = params.sv_path

    dmp_dir = params.dmp_dir
    mirror_bam_dir = params.mirror_bam_dir
    mirror_access_bam_dir = params.mirror_access_bam_dir
    dmp_key_path = file(params.dmp_key_path)
    access_key_path = file(params.access_key_path)

    GenerateBamPaths(
        samples = samples_ch, 
        id_mapping = id_ch, 
        bam_path_normal = bam_path_normal,
        bam_path_plasma_duplex = bam_path_plasma_duplex,
        bam_path_plasma_simplex = bam_path_plasma_simplex,
        bam_path_index_normal = bam_path_index_normal,
        bam_path_index_duplex = bam_path_index_duplex,
        bam_path_index_simplex = bam_path_index_simplex,
        maf_path = maf_path,
        cna_path = cna_path,
        sv_path = sv_path,
        dmp_dir = dmp_dir,
        mirror_bam_dir = mirror_bam_dir,
        mirror_access_bam_dir = mirror_access_bam_dir,
        dmp_key_path = dmp_key_path,
        access_key_path = access_key_path)

}
