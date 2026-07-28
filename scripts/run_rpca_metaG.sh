#to run RPCA you need to be in gemelli-qiime2-amplicon-2024.2

#AD all
path=/mnt/zarrinpar/scratch/sfloresr/AD_EM/metagenomics/
mtb_file=$path/wol2_genome_none.gemelli.distance.qza
md_file=$path/overly_specified_metadata.tsv
output_path=$path/rpca_results_all/
comparison=genotype_diet
comparison2=zt_collection_time
comparison3=sex
#comparison4=
#comparison5=

#biom convert \
# -i $mtb_file \
# -o ${mtb_file/%.tsv/.biom} \
# --table-type="OTU table" \
# --to-hdf5

#qiime tools import \
# --input-path ${mtb_file/%.tsv/.biom} \
# --output-path ${mtb_file/%.tsv/.qza} \
# --type FeatureTable[Frequency]

#qiime gemelli rpca \
#  --i-table ${mtb_file/%.tsv/.qza} \
#  --p-min-feature-count 0 \
#  --p-min-sample-count 0 \
#  --o-biplot $output_path/ordination.qza \
#  --o-distance-matrix $output_path/distance_matrix.qza

qiime emperor biplot \
  --i-biplot wol2_genome_none.gemelli.biplot.qza \
  --m-sample-metadata-file $md_file \
  --p-ignore-missing-samples \
  --o-visualization $output_path/biplot.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $mtb_file \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $mtb_file \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison2 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison2}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $mtb_file \
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
