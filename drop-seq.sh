#!/bin/sh
#$ -S /bin/sh
#$ -l s_vmem=80G
#$ -l mem_req=80G
#-pe def_slot 4

<< COMMENT
USAGE:  qsub -cwd -l os7 ~/Dropseq_dir/drop-seq.sh \
	-1 ~/Dropseq_dir/data/R1_barcode.fastq.gz \
	-2 ~/Dropseq_dir/data/R2_transcript.fastq.gz \
	-i ID_NAME -r REFERENCE -c EC

    Desctiption:
	mapping singlecell data and making umi count data

    Options
        -1 ... Input Read1 path (barcode read)
        -2 ... Input Read2 path (transcript read)
        -i ... analysis ID. This ID is used output file names
	-r ... reference ID (hg38 or mm10)
        -c ... Expected detect Cell Num. Default is DROP = 1000
COMMENT

#modules(python, bowtie, star, R, samtools)
module use /usr/local/package/modulefiles
module load python/2.7
module load bowtie/2.3.4.3
module load star/2.6.1c
module load r/3.5
module load samtools/1.9

#tools(umi_tools, cellranger, homer, STAR, subread, bedtools2)
export PATH=~/.local/bin/:$PATH
export PATH=~/Dropseq_dir/tools/cellranger-3.0.2/:$PATH
export PATH=~/Dropseq_dir/tools/homer/bin/:$PATH
export PATH=~/Dropseq_dir/tools/STAR-2.7.1a/bin/Linux_x86_64/:$PATH
export PATH=~/Dropseq_dir/tools/subread-1.6.4-Linux-x86_64/bin/:$PATH
export PATH=~/Dropseq_dir/tools/bedtools2/:$PATH

#barcode pattern
DROP="CCCCCCCCCCCCNNNNNNNNN"
L_LIST=("DROP")
DEF_NUM=(1000)
NUM_THREAD=8

while getopts 1:2:i:r:c::zh OPT;do
        case $OPT in
                1 )R1=$OPTARG #read1 file
                        ;;
                2 )R2=$OPTARG #read2 file
                        ;;
                i )ID=$OPTARG #analysis ID
                        ;;
		r )REF=$OPTARG #reference
			;;
                c )CELL=$OPTARG #CELL NUM
                        ;;
                z )RETRY=1;
                        ;;
                h )usage
                        ;;
                *?)usage
                        ;;
        esac
done

#Read1 Read2
if [ -z "$R2" ] || [ -z "$R1" ] ;then
        echo "Read1かRead2が未定義です";
fi

if [ -e $R1 ] && [ -e $R2 ] ;then
        :
else
        echo -e "Read1かRead2のどちらかが存在していません";
        echo -e "Read1 -f: $R1";
        echo -e "Read2 -r: $R2";
        exit
fi

#reference
if [ $REF = "GRCh38" || $REF = "hg38" ];then
	RIB=~/Dropseq_dir/reference/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/humRibosomal
	STAR_REF=~/Dropseq_dir/reference/refdata-cellranger-GRCh38-3.0.0/star/
	GTF=~/Dropseq_dir/reference/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
elif [ $REF = "mm10" ]; then
	RIB=~/Dropseq_dir/reference/Mus_musculus/UCSC/mm10/Sequence/AbundantSequences/musRibosomal
	STAR_REF=~/Dropseq_dir/reference/refdata-cellranger-mm10-1.2.0/star/
	GTF=~/Dropseq_dir/reference/refdata-cellranger-mm10-1.2.0/genes/genes.gtf
else
	echo "specified reference $REF is not available"
fi

# ID CEHCK
tmp_id=${R1/.*/}
ID=${ID:-$tmp_id}

#GET_BARCODE
BC=`eval echo '$'$MODE`

#CELL NUM
if [ -z $CELL ] ; then
case "$MODE" in
        ${L_LIST[0]} )CELL=${DEF_NUM[0]} ;; #DROP
esac
fi

WDIR=work_${ID}

if [ -e "$WDIR/${ID}_START_LOG" ] && [ -z "$RETRY" ]; then
        printf "\tsame Project \"%s\"  already exists.\n\tplease change ID\n" $ID;
        exit 1;
fi

if [ ! -e $WDIR ]; then
    mkdir -p $WDIR
else
        if [ -z "$RETRY" ]; then
                printf "\tsame Project \"%s\"  already exists.\n\tplease change ID\n" $WDIR;
                exit 1;
        fi
fi

printf "STARTING SINGLE CELL PIPELINE...\n INPUT READ1 $R1 \n INPUT READ2 $R2 \n Library Type is $MODE \n analysis ID is $ID \n expected cell num is $CELL \n Barcode pattern is $BC \n" >> "$WDIR/${ID}_START_LOG"

pR1=`basename $R1`
pR2=`basename $R2`

barcode(){
        D1=$WDIR/whitelist
        mkdir $D1
        umi_tools whitelist --stdin $R1 \
                            --bc-pattern=$BC \
                            --set-cell-number=$CELL \
                            --plot-prefix=${D1}/${ID}_expect_whitelist \
                            --log2stderr > ${D1}/${ID}_whitelist.txt;
        printf "make whitelist finish...\n" >> "$WDIR/${ID}_START_LOG"

        D2=$WDIR/extract
        mkdir $D2
        umi_tools extract --stdin $R1 \
                          --bc-pattern=$BC \
                          --read2-in=$R2 \
                          --stdout=${D2}/processed_${pR1} \
                          --read2-out=${D2}/processed_${pR2} \
                          --filter-cell-barcode \
                          --whitelist=${D1}/${ID}_whitelist.txt
        printf "extract barcode finish...\n" >>"$WDIR/${ID}_START_LOG"
        return 0;
}

mapping(){
        D3=$WDIR/rmRIB
        mkdir $D3
        bowtie2 -p $NUM_THREAD --un-gz ${D3}/processed_rmRIB_${pR2}  -x $RIB -U ${D2}/processed_${pR2}  \
        -S ${WDIR}/${ID}_RIB.sam > ${D3}/${ID}_bowtie2_log 2>&1
        rm ${WDIR}/${ID}_RIB.sam &

        D4=$WDIR/STAR
        mkdir $D4
        STAR --genomeDir $STAR_REF \
             --readFilesIn ${D3}/processed_rmRIB_${pR2} \
             --readFilesCommand zcat \
             --outFilterMultimapNmax 1 \
             --outSAMtype BAM SortedByCoordinate \
             --outWigType wiggle read1_5p \
             --sjdbGTFfile $GTF \
             --runThreadN $NUM_THREAD \
             --outFileNamePrefix ${WDIR}/${ID}_

        mv ${WDIR}/${ID}_Log* ${WDIR}/${ID}_SJ.out.tab ${WDIR}/${ID}__STARgenome $D4
        printf "mapping finish...\n" >>"$WDIR/${ID}_START_LOG"
}

fcount(){
        D5=$WDIR/featureCounts
        mkdir $D5
        featureCounts -a $GTF \
                      -o ${D5}/gene_assigned \
                      -g gene_name \
                      -R BAM ${WDIR}/${ID}_Aligned.sortedByCoord.out.bam \
                      -s 1 \
                      -T 4

        samtools sort ${D5}/${ID}_Aligned.sortedByCoord.out.bam.featureCounts.bam ${WDIR}/${ID}_sorted
        samtools index ${WDIR}/${ID}_sorted.bam
}

umi_c(){
        D6=$WDIR/umi_count
        mkdir $D6
        umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --per-cell -I ${WDIR}/${ID}_sorted.bam -S ${WDIR}/${ID}_counts.tsv.gz

        mv ${WDIR}/${ID}_sorted.bam* ${WDIR}/umi_count
        printf "umi count finish...\n" >>"$WDIR/${ID}_START_LOG"
}

bigwig(){
        #D7=$WDIR/homer
        #mkdir $D7
        makeTagDirectory ${WDIR}/bigwig $D6/${ID}_sorted.bam -keepAll -mapq -1
        makeUCSCfile ${WDIR}/bigwig/ -o auto -noadj
        gunzip ${WDIR}/bigwig/bigwig.ucsc.bedGraph.gz
        bedGraphToBigWig ${WDIR}/bigwig/bigwig.ucsc.bedGraph $hg38 ${WDIR}/bigwig/bigwig.ucsc.bw
	printf "transformation to bigwig finish...\n" >>"$WDIR/${ID}_START_LOG"
}

#MAIN
barcode
mapping
fcount
umi_c
bigwig
exit;
