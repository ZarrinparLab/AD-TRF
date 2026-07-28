#to run RPCA you need to be in gemelli-qiime2-amplicon-2024.2

#AD all
path=/mnt/zarrinpar/scratch/sfloresr/AD_EM/metabolomics/
mtb_file=$path/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_feces.txt
md_file=$path/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-forrpca.txt
output_path=$path/diversity_metrics/rpca_fecal_all/
comparison=ATTRIBUTE_geno_n_feedingcond
comparison2=genofedZT
comparison3=ATTRIBUTE_ZT_collection_time
comparison4=BiologicalSex
#comparison5=

biom convert \
-i $mtb_file \
-o ${mtb_file/%.txt/.biom} \
--table-type="OTU table" \
--to-hdf5

qiime tools import \
--input-path ${mtb_file/%.txt/.biom} \
--output-path ${mtb_file/%.txt/.qza} \
--type FeatureTable[Frequency]

qiime gemelli rpca \
 --i-table ${mtb_file/%.txt/.qza} \
 --p-min-feature-count 0 \
 --p-min-sample-count 0 \
 --o-biplot $output_path/ordination.qza \
 --o-distance-matrix $output_path/distance_matrix.qza

qiime emperor biplot \
  --i-biplot $output_path/ordination.qza \
  --m-sample-metadata-file $md_file \
  --p-ignore-missing-samples \
  --o-visualization $output_path/biplot.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison2 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison2}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison3 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison3}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison4 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison4}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison5 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison5}-significance.qzv
