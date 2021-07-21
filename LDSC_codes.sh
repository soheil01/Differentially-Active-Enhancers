
https://github.com/bulik/ldsc

## Active LDSC environment
source activate ldsc

## All pCR coordinates
pCR_202163_sorted.bed

#Step:1
### Make annotation file
for chr in {1..22}

do

python ldsc/make_annot.py \
--bed-file LDSC_analysis/pCR_202163_sorted.bed \
--bimfile LDSC_analysis/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
--annot-file LDSC_analysis/Annotation.LD_files2/pCR/pCR_202163.${chr}.annot.gz

done

# Step 2:
### calculate LDSC
for chr in {1..22}

do

python ldsc/ldsc.py \
--l2 \
--bfile LDSC_analysis/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
--ld-wind-cm 1 \
--thin-annot \
--annot LDSC_analysis/Annotation.LD_files2/pCR/pCR_202163.${chr}.annot.gz \
--out LDSC_analysis/Annotation.LD_files2/pCR/pCR_202163.${chr} \
--print-snps LDSC_analysis/hapmap3_snps/ 

done

#Step3:
### Convert GWAS summary to proper .sumstats format
python LDSC_analysis/ldsc/munge_sumstats.py \
--sumstats Autism_Grove_iPSYCH-PGC_ASD_Nov2017.gz \  ## replace based on GWAS data
--out LDSC_analysis/Reformatting_Summary_Statistics_GWASData/Autism_Grove_2019 \
--merge-alleles LDSC_analysis/w_hm3.snplist \

#Step4:
### calculate partitaion heritability
python LDSC_analysis/ldsc/ldsc.py \
--h2 LDSC_analysis/Reformatting_Summary_Statistics_GWASData/Autism_Grove_2019.sumstats.gz \  ## replace based on GWAS data
--ref-ld-chr LDSC_analysis/Annotation.LD_files2/pCR/pCR_202163.,LDSC_analysis/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr LDSC_analysis/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr LDSC_analysis/1000G_Phase3_frq/1000G.EUR.QC. \
--print-coefficients \
--out LDSC_analysis/partitaion_heritability_results/pCR/Autism_Grove_2019_baseline



