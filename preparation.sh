#!/bin/sh
#$ -S /bin/sh
#$ -l s_vmem=48G
#$ -l mem_req=48G
#-pe def_slot 4

<<COMMENT
usage: qsub -l os7 -cwd preparation.sh
execute in your Dropseq directory
COMMENT

mkdir tools
cd tools

wget -O cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1562278360&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU2MjI3ODM2MH19fV19&Signature=HW0lNZbObneP5KRm-5CAsoRkvU2IjVNI4r-5vxtNDtxqmnFJvsNFw6it6PCP16vIFF89N73Kh-3Ss6djZUy2vxE4RhEkER0PgYxzbUxN5AW~Qr9HLlSEqnpY5ATnF-nP4mJrx7mxZvc2tY4Dm3tJdGg1uyjdnk7K6565tf-R2vaEhcKBB55jdvF9PIjBWtMGSNkZFzXdd-8jpT2PXD1dtunEZKYyveBcBcb6dFeBvQ-bmcbO3PFjUgj5Un0M6vJtEkvPS3vmWmmGm1rUfBgrHf5oOtIF7nOT2BXZiYbSta1cefKhxZ6mVXoEPZf07gvyaBr7oUeT8QMb5frKmHSIyg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -xzvf cellranger-3.0.2.tar.gz

mkdir homer
cd homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
cd ../

wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
tar -zxvf bedtools-2.28.0.tar.gz
cd bedtools2
make
cd ../

wget https://github.com/alexdobin/STAR/archive/2.7.1a.tar.gz
tar -xzf 2.7.1a.tar.gz
cd STAR-2.7.1a/source
make
cd ../../

wget https://sourceforge.net/projects/subread/files/subread-1.6.4/subread-1.6.4-Linux-x86_64.tar.gz
tar zxvf subread-1.6.4-Linux-x86_64.tar.gz
cd ../

export PATH=~/Dropseq_dir/tools/cellranger-3.0.2:$PATH
mkdir reference
cd reference

wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz

wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
gunzip Homo_sapiens.GRCh38.93.gtf.gz


cellranger mkgtf Homo_sapiens.GRCh38.93.gtf Homo_sapiens.GRCh38.93.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene


cellranger mkref --genome=GRCh38 \
                 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh38.93.filtered.gtf \
                 --ref-version=3.0.0

wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-3.0.0.tar.gz

wget ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz


wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
gunzip Mus_musculus.GRCm38.93.gtf.gz


cellranger mkgtf Mus_musculus.GRCm38.93.gtf Mus_musculus.GRCm38.93.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene


cellranger mkref --genome=mm10 \
                 --fasta=Mus_musculus.GRCm38.dna.primary_assembly.fa \
                 --genes=Mus_musculus.GRCm38.93.filtered.gtf \
                 --ref-version=3.0.0

wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
tar xvf Homo_sapiens_UCSC_hg38.tar.gz

wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar xvf Mus_musculus_UCSC_mm10.tar.gz

cd ../
pip install --user umi_tools
















