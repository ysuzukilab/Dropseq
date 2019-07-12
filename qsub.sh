#!/bin/sh
#$ -S /bin/sh

<< COMMENT
usage:
	sh qsub.sh -r REFERENCE -f FASTQS -e EXPECTCELLS
description:
	If you have more than one fastq dataset in data directory and 
	want to run all of them at the same time, write the path to
	data directory after -f. Otherwise, write the path to fastqs.
options:
        -r ... reference (GRCh38/hg38 or mm10)
	-f ... path to data directory or fastqs
        -e ... expected number of cells
COMMENT

while getopts r:f:e: OPT;do
        case $OPT in
                r )REF=$OPTARG #REFERENCE
                        ;;
                f )FASTQS=$OPTARG #path to data directory or fastqs
                        ;;
                e )EC=$OPTARG #expected number of cells
                        ;;
		*?)usage
                        ;;
        esac
done

if [ -z "$REF" ];then
        echo "Undefined REFERENCE";
        exit 1;
fi

if [ -z "$FASTQS" ];then
        echo "path to FASTQS does not exist";
        exit 1;
fi

if [ -z "$EC" ];then
        echo "Undefined EXPECTCELLS";
        exit 1;
fi

R1S=($(ls -1 ${FASTQS}/*R1*))
R2S=($(ls -1 ${FASTQS}/*R2*))
NUM=($(ls -1 ${FASTQS}/*R1* | wc -l))

for i in `seq 0 $(($NUM - 1))`
do
	qsub -l os7 -cwd ~/Dropseq_dir/drop-seq.sh -1 ${R1S[$i]} -2 ${R2S[$i]} -i sample$(($i + 1)) -r $REF -c $EC
done

