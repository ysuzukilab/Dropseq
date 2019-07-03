#!/bin/sh
#$ -S /bin/sh

<< COMMENT
USAGE: qsub -cwd -l s_vmem=50G,mem_req=50G ~/Dropseq_dir/drop-seq.sh -f ~/Dropseq_dir/data/R1_barcode.fastq.gz -r ~/Dropseq_dir/data/R2_transcript.fastq.gz -i ID_NAME -c EC -m DROP

    Desctiption:
          mapping singlecell data and making umi count data

    Options
        -f ... Input Read1 path (barcode read)
        -r ... Input Read2 path (transcript read)
        -i ... analysis ID .This ID use output file names
        -c ... Expected detect Cell Num. default is
               DROP      = 500
        -m ... The following Library is currently available.
               "DDSEQ" "DROP"
COMMENT

#tools(python, bowtie2, star, subread, R)
module use /usr/local/package/modulefiles
module load 

#refence(mouse)

#refence(human)

#barcode pattern
DDSEQ=".{0,5}(?P<cell_1>.{6})(?P<discard_1>TAGCCATCGCATTGC)(?P<cell_2>.{6})(?P<discard_2>TACCTCTGAGCTGAA)(?P<cell_3>.{6})(?P<discard_2>ACG)(?P<umi_1>.{8})GAC.*"
ICELL8=""
DROP="CCCCCCCCCCCCNNNNNNNNN"

L_LIST=("DROP")
DEF_NUM=(1000)

NUM_THREAD=8

while getopts f:r:i:c:m::zh OPT;do
        case $OPT in
                f )R1=$OPTARG #read1 file
                        ;;
                r )R2=$OPTARG #read1 file
                        ;;
                i )ID=$OPTARG #analysis ID
                        ;;
                c )CELL=$OPTARG #CELL NUM
                        ;;
                m )MODE=$OPTARG #LIBRARY TYPE
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
        #usage
fi

if [ -e $R1 ] && [ -e $R2 ] ;then
        :
else
        echo -e "Read1かRead2のどちらかが存在していません";
        echo -e "Read1 -f: $R1";
        echo -e "Read1 -r: $R2";
        exit
fi

#mode
if ! `echo ${L_LIST[@]} | grep -wq "$MODE"` ; then
        echo "$COMM";
        exit 1;
fi


# ID CEHCK
tmp_id=${R1/.*/}
ID=${ID:-$tmp_id}

#GET_BARCODE
BC=`eval echo '$'$MODE`

#CELL NUM
if [ -z $CELL ] ; then
case "$MODE" in
        ${L_LIST[0]} )CELL=${DEF_NUM[0]} ;; #DDSEQ
        #${L_LIST[1]} )CELL=${DEF_NUM[1]} ;; #ICELL8
        #${L_LIST[2]} )CELL=${DEF_NUM[2]} ;; #DROP
        #${L_LIST[3]} )CELL=${DEF_NUM[3]} ;; #ICELL8
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

<<Commentout
Rscript(){
        R_LIBS_USER=/sshare4/home/tissei/R/x86_64-pc-linux-gnu-library/3.5 /usr/local/package/r/3.5.0/bin/Rscript /home/tissei/rscript_graphs_drop.R ${WDIR}/${ID}_counts.tsv.gz ${WDIR}/analysis_summary.txt ${WDIR}/rscript.rds > ${WDIR}/analysis_summary.txt
	mv *.png ${WDIR}
	printf "statistical analysis with R finish...\n" >>"$WDIR/${ID}_START_LOG"
}
Commentout

#MAIN
barcode
mapping
fcount
umi_c
bigwig
#Rscript
exit;
