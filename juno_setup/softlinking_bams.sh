#!/usr/bin/env bash

export test_bam_dir=/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/bams
rm -rf $test_bam_dir/*

ln -s /juno/work/access/production/runs/Project_06302_?/bam_qc/current/*/*/C-AKJWE7*.ba? $test_bam_dir
ln -s /juno/work/access/production/runs/Project_06302_W/bam_qc/06302_W/*/*/C-AKJWE7*ba? $test_bam_dir
ln -s /juno/work/access/production/runs/Project_06302_?/bam_qc/current/*/*/C-X8EA4H*.ba? $test_bam_dir