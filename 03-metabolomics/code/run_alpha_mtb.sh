#to run RPCA you need to be in gemelli-qiime2-amplicon-2024.2

#AD all
path=/mnt/zarrinpar/scratch/sfloresr/AD_EM/metabolomics/
mtb_file=$path/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_feces.txt
md_file=$path/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-forrpca.txt
alpha_type=shannon 
output_path=$path/diversity_metrics/alpha_shannon_fecal_all/

biom convert \
-i $mtb_file \
-o ${mtb_file/%.txt/.biom} \
--table-type="OTU table" \
--to-hdf5

qiime tools import \
--input-path ${mtb_file/%.txt/.biom} \
--output-path ${mtb_file/%.txt/.qza} \
--type FeatureTable[Frequency]

qiime diversity alpha \
--i-table ${mtb_file/%.txt/.qza} \
--p-metric $alpha_type \
--o-alpha-diversity $output_path/alpha.qza

qiime tools export \
--input-path $output_path/alpha.qza \
--output-path $output_path/alpha