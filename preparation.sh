#!/bin/sh
#$ -S /bin/sh
#$ -l s_vmem=48G
#$ -l mem_req=48G
#-pe def_slot 4

<<COMMENT
usage: qsub -l os7 -cwd preparation.sh
execute in your 10X directory
COMMENT

mkdir tools
cd tools
wget -O cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1562030772&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU2MjAzMDc3Mn19fV19&Signature=eotQfzfJo6ebf6Y~erxL9qxtpoMDiQn2tw3r4RRT2JspqS4lOuST9HCAf1OGAOfiVWxLuUxPXJJVzztthX0Wwh7k2~MNvhyGNBkLgkERXUjvtAP6PsQEzM8xV~aGbpllbI63Dz0vgA~hLArwVg62jKhU3Jiz-KEtxu3jxvJtgjIoAL6wQa-ZGb4Z-fB~YXDGRcEkuAlQJy70cHl0qK42RN0ARsIiTvY~-nMNYGucVof~4WD3RieBoV2kWPe78mGPqk0DQNOmvqwp18NlXrCcmMnJM2CiJgv3rUyELmEdkN-zvdwHbFZasVOi-7iZymVa6T-I7iUOWurLzrAZjiANUQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
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




