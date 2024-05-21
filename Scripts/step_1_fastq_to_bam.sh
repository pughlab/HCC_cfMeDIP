# Description: Processes raw FASTQ files for genomic analysis, including FASTQ quality control, alignment, sorting, indexing, deduplication, and QC of deduplicated BAM files.
# Author：Kui Chen → kui.chen@uhn.ca
# 1st Created: 2021-10-01
# Last Modified: 2023-12-10

echo "Processing of raw FASTQ initiated......"

# Navigate to the directory containing the FASTQ files
cd /path/to/your/fastq/files

# Define directories
FASTQC_PRE_DIR="./fastqc_pre"
FASTQC_PST_DIR="./fastqc_pst"
FASTP_DIR="./fastp_report"
BAM_DIR="../bam"
QC_DIR="../fastq_bam_QC"
HG19_DIR="./hg19"

# Define and create necessary directories
mkdir -p $FASTQC_PRE_DIR $FASTQC_PST_DIR $FASTP_DIR $BAM_DIR $QC_DIR $BAM_QC_DIR $HG19_DIR

# Create soft-links for the hg19 reference files in the hg19 directory
ln -s /path/to/original/hg19.1.bt2 $HG19_DIR/hg19.1.bt2
ln -s /path/to/original/hg19.2.bt2 $HG19_DIR/hg19.2.bt2
ln -s /path/to/original/hg19.3.bt2 $HG19_DIR/hg19.3.bt2
ln -s /path/to/original/hg19.4.bt2 $HG19_DIR/hg19.4.bt2
ln -s /path/to/original/hg19.rev.1.bt2 $HG19_DIR/hg19.rev.1.bt2
ln -s /path/to/original/hg19.rev.2.bt2 $HG19_DIR/hg19.rev.2.bt2
ln -s /path/to/original/hg19.fa $HG19_DIR/hg19.fa
ln -s /path/to/original/hg19.fa.fai $HG19_DIR/hg19.fa.fai
ln -s /path/to/original/hg19.gtf $HG19_DIR/hg19.gtf

# Reference genome path
REF_INDEX="$HG19_DIR/hg19"
GTF_FILE="$HG19_DIR/hg19.gtf"

# Merge fastq files from libraries split across different lanes
for file in *_L001_R1_001.fastq.gz; do
    base=$(basename $file _L001_R1_001.fastq.gz)
    if [[ -f "${base}_L002_R1_001.fastq.gz" ]]; then
        # Merge R1 and R2 reads
        cat ${base}_L001_R1_001.fastq.gz ${base}_L002_R1_001.fastq.gz > ${base}_R1_001.fastq.gz
        cat ${base}_L001_R2_001.fastq.gz ${base}_L002_R2_001.fastq.gz > ${base}_R2_001.fastq.gz

        # Remove original files if merge is successful
        rm ${base}_L001_R1_001.fastq.gz ${base}_L002_R1_001.fastq.gz ${base}_L001_R2_001.fastq.gz ${base}_L002_R2_001.fastq.gz
    fi
done

echo "Fastq files merging completed."

# Perform FastQC on all R1 and R2 files after merging
for file in *_R1_001.fastq.gz; do
    base=$(basename $file _R1_001.fastq.gz)
    fastqc -t 10 -o $FASTQC_PRE_DIR ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz
done

echo "FastQC preparation completed."

# Main processing loop for each sample
for R1 in *_R1_001.fastq.gz; do
    SAMPLE=$(basename $R1 _R1_001.fastq.gz)
    R2=${SAMPLE}_R2_001.fastq.gz

    # fastp performs adapter trimming, quality filtering, and generates HTML reports for quality control
    fastp -w 10 -i $R1 -I $R2 -o ${SAMPLE}_fp_R1.fastq.gz -O ${SAMPLE}_fp_R2.fastq.gz -g -f 5 -F 5 --correction --html=${SAMPLE}_fastp.html
    mv ${SAMPLE}_fastp.html $FASTP_DIR
    rm *.json
    # IFSD (Insert Freagments Size Distribution) table were extract form the fastp html reoprt to generate Figure S2B

    # Perform post-processed QC with FastQC
    fastqc -t 10 -o $FASTQC_PST_DIR ${SAMPLE}_fp_R1.fastq.gz ${SAMPLE}_fp_R2.fastq.gz

    # Align reads to the reference genome with Bowtie2
    bowtie2 -p 10 --minins 50 --maxins 500 -x $REF_INDEX -1 ${SAMPLE}_fp_R1.fastq.gz -2 ${SAMPLE}_fp_R2.fastq.gz -S ${SAMPLE}.sam

    # Convert SAM to BAM, sort, and index
    samtools view -@ 10 -bS -f 2 -F 4 ${SAMPLE}.sam > ${SAMPLE}.bam
    samtools sort -o ${SAMPLE}.sort.bam -@ 10 ${SAMPLE}.bam
    samtools index ${SAMPLE}.sort.bam

    # Rename sorted BAM files for consistency
    mv ${SAMPLE}.sort.bam ${SAMPLE}.bam
    mv ${SAMPLE}.sort.bam.bai ${SAMPLE}.bam.bai

    # Perform pre-QC for sorted BAM files
    samtools flagstat ${SAMPLE}.bam > ${SAMPLE}_pre.bam.flagstat.txt
    qualimap bamqc -bam ${SAMPLE}.bam -c --java-mem-size=30G -gff $GTF_FILE -outfile ${SAMPLE}_preQC -outformat PDF -nt 10

    # Deduplication using sambamba
    sambamba markdup -r -p -t 10 ${SAMPLE}.bam ${BAM_DIR}/${SAMPLE}_dedup.bam

    # Post-deduplication QC
    samtools flagstat ${BAM_DIR}/${SAMPLE}_dedup.bam > ${BAM_QC_DIR}/${SAMPLE}_dedup.bam.flagstat.txt
    qualimap bamqc -bam ${BAM_DIR}/${SAMPLE}_dedup.bam -c --java-mem-size=30G -gff $GTF_FILE -outfile ${BAM_QC_DIR}/${SAMPLE}_dedupQC -outformat PDF -nt 10

    # Rename deduplicated BAM files for simplicity
    mv ${BAM_DIR}/${SAMPLE}_dedup.bam ${BAM_DIR}/${SAMPLE}.bam
    mv ${BAM_DIR}/${SAMPLE}_dedup.bam.bai ${BAM_DIR}/${SAMPLE}.bam.bai

    # Merge all QC result
    mv ${FASTQC_PRE_DIR} ${QC_DIR}/fastqc_pre
    mv ${FASTQC_PST_DIR} ${QC_DIR}/fastqc_post
    mv ${FASTP_DIR} ${QC_DIR}/fastp_report
    mv ${BAM_QC_DIR} ${QC_DIR}/bam_qc

    # Delete intermediate fastp output files as they are no longer needed
    rm ${SAMPLE}_fp_R1.fastq.gz ${SAMPLE}_fp_R2.fastq.gz

    # Delete the intermediate SAM, BAM and its index file as they are no longer needed
    rm ${SAMPLE}.sam ${SAMPLE}.bam ${SAMPLE}.bam.bai

    echo "${SAMPLE} processing completed."

done

echo "All BAM files are ready!"
