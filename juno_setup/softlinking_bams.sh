#!/usr/bin/env bash

export test_dir=/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/
test_bam_dir="${test_dir}bams"
test_maf_dir="${test_dir}mafs"
test_cna_dir="${test_dir}cnas"
test_sv_dir="${test_dir}svs"

mkdir $test_dir
mkdir $test_bam_dir
mkdir $test_maf_dir
mkdir $test_cna_dir
mkdir $test_sv_dir
rm -rf $test_bam_dir/*
rm -rf $test_maf_dir/*
rm -rf $test_cna_dir/*
rm -rf $test_sv_dir/*
ls $test_dir/*

ln -s /juno/work/access/production/runs/Project_06302_?/bam_qc/current/*/*/C-AKJWE7*.ba? $test_bam_dir
ln -s /juno/work/access/production/runs/Project_06302_W/bam_qc/06302_W/*/*/C-AKJWE7*ba? $test_bam_dir
ln -s /juno/work/access/production/runs/Project_06302_?/bam_qc/current/*/*/C-X8EA4H*.ba? $test_bam_dir

cp /juno/work/bergerm1/bergerlab/zhengy1/TDM1/snv_pipeline/maf_101019/C-AKJWE7* $test_maf_dir
cp /juno/work/bergerm1/bergerlab/zhengy1/TDM1/snv_pipeline/maf_101019/C-X8EA4H* $test_maf_dir

head -1 /juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/cna_032120/TDM1_copynumber_segclusp.genes.txt > "${test_cna_dir}/cna_results.txt"
grep -P 'AKJWE7|X8EA4H' /juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/cna_032120/TDM1_copynumber_segclusp.genes.txt >> "${test_cna_dir}/cna_results.txt"

cp /juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/manta_032220/final_results/C-AKJWE7* $test_sv_dir
cp /juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/manta_032220/final_results/C-X8EA4H* $test_sv_dir
