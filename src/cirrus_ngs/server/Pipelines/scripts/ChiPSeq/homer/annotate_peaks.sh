#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
chip_sample=$5
input_sample=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped

#logging
log_dir=$log_dir/$chip_sample
mkdir -p $log_dir
log_file=$log_dir/'annotate_peaks.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$chip_sample
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

rm -r $workspace/tags_$chip_sample $workspace/tags_$input_sample &>/dev/null

pair_base_name=$chip_sample

##DOWNLOAD##
if [ ! -f $workspace/$chip_sample$style_ext ] || [ "$pairs_exist" == "True" ]
then
    if [ "$pairs_exist" == "True" ]
    then
        pair_base_name=$chip_sample'_vs_'$input_sample
        input_address=$input_address/$chip_sample
    fi

    aws s3 cp $input_address/$pair_base_name$style_ext $workspace/ --quiet
fi
##END_DOWNLOAD##


##ANNOTATEPEAKS##
mkdir -p $workspace/go_$pair_base_name
mkdir -p $workspace/ontology_$pair_base_name

check_exit_status "$annotate_peaks $workspace/$pair_base_name$style_ext $genome \
    -go "$workspace/go_$pair_base_name" -genomeOntology "$workspace/ontology_$pair_base_name" \
    > $workspace/$pair_base_name.annotated$style_ext" $JOB_NAME $status_file

go_dir_outputs=($workspace/go_$pair_base_name/biocyc.txt
$workspace/go_$pair_base_name/biological_process.txt
$workspace/go_$pair_base_name/cellular_component.txt
$workspace/go_$pair_base_name/chromosome.txt
$workspace/go_$pair_base_name/cosmic.txt
$workspace/go_$pair_base_name/gene3d.txt
$workspace/go_$pair_base_name/geneOntology.html
$workspace/go_$pair_base_name/gwas.txt
$workspace/go_$pair_base_name/interactions.txt
$workspace/go_$pair_base_name/interpro.txt
$workspace/go_$pair_base_name/kegg.txt
$workspace/go_$pair_base_name/lipidmaps.txt
$workspace/go_$pair_base_name/molecular_function.txt
$workspace/go_$pair_base_name/msigdb.txt
$workspace/go_$pair_base_name/pathwayInteractionDB.txt
$workspace/go_$pair_base_name/pfam.txt
$workspace/go_$pair_base_name/prints.txt
$workspace/go_$pair_base_name/prosite.txt
$workspace/go_$pair_base_name/reactome.txt
$workspace/go_$pair_base_name/smart.txt
$workspace/go_$pair_base_name/smpdb.txt
$workspace/go_$pair_base_name/wikipathways.txt)

ontology_dir_outputs=($workspace/ontology_$pair_base_name/basic.genomeOntology.txt
$workspace/ontology_$pair_base_name/GenomeOntology.html
$workspace/ontology_$pair_base_name/repeats.genomeOntology.txt)

check_exit_status "check_outputs_exist ${go_dir_outputs[@]} ${ontology_dir_outputs[@]} $workspace/$pair_base_name.annotated$style_ext" $JOB_NAME $status_file
##END_ANNOTATEPEAKS##


##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$pair_base_name.annotated$style_ext" --recursive --quiet
aws s3 cp $workspace/go_$pair_base_name $output_address/go_$pair_base_name --recursive --quiet
aws s3 cp $workspace/ontology_$pair_base_name $output_address/ontology_$pair_base_name --recursive --quiet
##END_UPLOAD##
