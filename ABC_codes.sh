## Active ABC environment
conda activate final-abc-env

## The below file prepared in R. Extending the regions 500 bps (250bps from each side)of center
PCR_202163_ABCInteractin_500bp_Included.CollapsedGeneBounds.TSS500bp.bed

K27ac_data = Data/K27ac/ # Including all relevant K27ac 
DNase_data = Data/DNase/  # Including all relevant DNase 

## Quantifying Enhancer Activity:
python src/run.neighborhoods.py \
--candidate_enhancer_regions ABC/Input/PCR_202163_ABCInteractin_500bp_Included.CollapsedGeneBounds.TSS500bp.bed \
--genes reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed \
--H3K27ac  K27ac_data
--DHS DNase_data
--expression_table ABC/Input/Fetal.Brain_Roadmap.txt \
--chrom_sizes reference/chr_sizes \
--ubiquitously_expressed_genes reference/UbiquitouslyExpressedGenesHG19.txt \
--outdir ABC/output/Neighborhoods/

## Predict target genes for CP and GZ HiC data
python src/predict.py \
--enhancers ABC/output/Neighborhoods/EnhancerList.txt \
--genes ABC/output/Neighborhoods/GeneList.txt \
--HiCdir ABC/Input/HiC_FetalBrain_Won/CP  \
--chrom_sizes reference/chr_sizes \
--hic_type bedpe \
--hic_resolution 10000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--outdir ABC/output/Predictions/CP_interactions.02/ \
--make_all_putative


python src/predict.py \
--enhancers ABC/output/Neighborhoods/EnhancerList.txt \
--genes ABC/output/Neighborhoods/GeneList.txt \
--HiCdir ABC/Input/HiC_FetalBrain_Won/GZ  \
--chrom_sizes reference/chr_sizes \
--hic_type bedpe \
--hic_resolution 10000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--outdir ABC/output/Predictions/GZ_interactions.02/ \
--make_all_putative
