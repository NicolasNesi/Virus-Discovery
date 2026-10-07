#!/bin/bash
#SBATCH -J DetectionVirusP1                     # Job name 
#SBATCH -o DetectionVirusP1.%a.%A.out           # File to which stdout will be written
#SBATCH -e DetectionVirusP1%a.%A.err            # File to which stderr will be written
#SBATCH -p long                                 # Partition
#SBATCH -N 1                                    # Number of node
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH --cpus-per-task=48                      # Number of cpu per task
#SBATCH -t 04-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 500G                              # Memory for all cores in Mbytes (--mem-per-cpu for MPI jobs)
#SBATCH --mail-type=ALL                         # BEGIN,END,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr     # Email address

# ---------------------------------
# Version v1.0
# ---------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 05/05/2025
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module purge
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate MinerVirus_env
# include:
# blast=2.14.1 
# bowtie2=2.4.5
# diamond=2.1.8
# fastp=0.23.4
# jellyfish=2.2.10 
# megahit=1.2.9
# rsem=1.3.3
# salmon=0.14.2
# samtools=1.18
# trinity=2.15.1
# ---------------------------------

# ---------------------------------
# Variables paths
# ---------------------------------
GROUP="${1}"
PROJECT="${2}"
FOLDER="/home/2019013/Data/$GROUP/$PROJECT"
inpath="${FOLDER}/raw_reads"
trimmed="${FOLDER}/trimmed_reads"
outpath="${FOLDER}/contigs"
outpath2="${FOLDER}/blast_results"
outfinal="${FOLDER}/contigs/final_contigs"
outabundance="${FOLDER}/abundance/final_abundance"
tempdir="${FOLDER}/blast_results"
outlogs="${FOLDER}/logs"

# Database for BLAST
db1="/home/2019013/Databases/RdRp-scan/RdRp-scan_0.90.dmnd"
db2="/home/2019013/Databases/RVDB/U-RVDBv26.0-prot.fasta"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" ${FOLDER}/Scripts/List_samples_${RUN}.txt)"
RUNDATE="$(date '+%Y-%m-%d')"

echo "Job started with SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
# Create directories
# ---------------------------------
# List of directories to create
dirs=("$outpath" "$trimmed" "$outfinal" "$outabundance" "$outlogs")

# Loop through and create if they don't exist
for dir in "${dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi
done
# ---------------------------------

# ---------------------------------
# Quality trimming Illumina
# ---------------------------------
echo "Starting reads trimming and filtering using Fastq for sample $SAMPLE"

fastp \
--in1 "${inpath}/${SAMPLE}_R1.fastq" \
--in2 "${inpath}/${SAMPLE}_R2.fastq" \
--detect_adapter_for_pe \
--qualified_quality_phred 15 \
--low_complexity_filter \
--dedup \
--thread 24 \
--json "${SAMPLE}"_ReportFastp.json \
--html "${SAMPLE}"_ReportFastp.html \
--report_title "${SAMPLE} Illumina Fastp report" \
--out1 "${trimmed}/${SAMPLE}_trimmed_R1.fastq" \
--out2 "${trimmed}/${SAMPLE}_trimmed_R2.fastq" 2>&1 | tee "${SAMPLE}_fastp.log"
# ---------------------------------

# ---------------------------------
# Assembly Megahit
# ---------------------------------
echo "Starting reads denovo assembly using megahit for sample $SAMPLE"

 megahit \
-1 "${trimmed}/${SAMPLE}_trimmed_R1.fastq" \
-2 "${trimmed}/${SAMPLE}_trimmed_R2.fastq" \
--num-cpu-threads "${THREADS}" \
--memory 0.9 \
-o "${outpath}/${SAMPLE}_out"

cat "${outpath}/${SAMPLE}_out/final.contigs.fa" | sed "s/=//g" | sed "s/ /_/g" > "${outpath}/${SAMPLE}.contigs.fa"
cp "${outpath}/${SAMPLE}_out/log" "${outlogs}/${SAMPLE}.log"
rm -r "${outpath}/${SAMPLE}_out"
# ---------------------------------

# ---------------------------------
# Abundance Trinity
# ---------------------------------
echo "Starting quantification of abundance using RSEM for $SAMPLE"

align_and_estimate_abundance.pl \
--transcripts "${outpath}/${SAMPLE}.contigs.fa" \
--seqType fq \
--left "${trimmed}/${SAMPLE}_trimmed_R1.fastq" \
--right "${trimmed}/${SAMPLE}_trimmed_R2.fastq" \
--est_method RSEM \
--aln_method bowtie2 \
--output_dir "${FOLDER}/abundance/${SAMPLE}_abundance" \
--thread_count "${THREADS}" \
--prep_reference

cp "${FOLDER}/abundance/${SAMPLE}_abundance/RSEM.isoforms.results" "${outabundance}/${SAMPLE}_RSEM.isoforms.results"
#rm $FOLDER/abundance/${SAMPLE}_abundance/bowtie2.bam
#rm $outfinal/${SAMPLE}.contigs.fa.*
# ---------------------------------

# ---------------------------------
# Diamond blastx versus RdRp databases
# ---------------------------------
echo "Starting Diamond Blastx against RdRp database for $SAMPLE"

diamond blastx \
--query "${outpath}/${SAMPLE}.contigs.fa" \
--db "${db1}" \
--tmpdir "${tempdir}" \
--out "${outpath2}/${SAMPLE}_rdrp_blastx_results.txt" \
--evalue 1E-4 \
-c2 \
--max-target-seqs 3 \
--block-size 14 \
--threads "${THREADS}" \
--outfmt 6 qseqid qlen sseqid stitle pident length evalue \
--ultra-sensitive
# ---------------------------------

# ---------------------------------
# Extract contigs from BLAST to fasta
# ---------------------------------
grep -i ".*" "${outpath2}/${SAMPLE}_rdrp_blastx_results.txt" | cut -f1 | sort | uniq > "${outpath2}/${SAMPLE}_temp_contig_names.txt" 
grep -A1 -I -Ff "${outpath2}/${SAMPLE}_temp_contig_names.txt" "${outpath}/${SAMPLE}.contigs.fa" > "${outpath2}/${SAMPLE}_rdrp_blastcontigs.fasta"
sed -i 's/--//' "${outpath2}/${SAMPLE}_rdrp_blastcontigs.fasta"
sed -i '/^[[:space:]]*$/d' "${outpath2}/${SAMPLE}_rdrp_blastcontigs.fasta"
sed --posix -i "/^\>/ s/$/"_$SAMPLE"/" "${outpath2}/${SAMPLE}_rdrp_blastcontigs.fasta"
rm "${outpath2}/${SAMPLE}_temp_contig_names.txt" 
# ---------------------------------

# ---------------------------------
# Diamond blastx versus RVRB databases
# ---------------------------------
echo "Starting Diamond Blastx against RVRB database for $SAMPLE"

diamond blastx \
--query "${outpath}/${SAMPLE}.contigs.fa" \
--db "${db2}" \
--tmpdir "${tempdir}" \
--out "${outpath2}/${SAMPLE}_RVRB_blastx_results.txt" \
--evalue 1E-10 \
-c1 \
--max-target-seqs 1 \
--block-size 14 \
--threads "${THREADS}" \
--outfmt 6 qseqid qlen sseqid stitle pident length evalue \
--ultra-sensitive \
--iterate
# ---------------------------------

# ---------------------------------
# Extract contigs
# ---------------------------------
grep -i ".*" "{$outpath2}/${SAMPLE}_RVRB_blastx_results.txt" | cut -f1 | sort | uniq > "${outpath2}/${SAMPLE}_temp_contig_names.txt" 
grep -A1 -I -Ff "${outpath2}/${SAMPLE}_temp_contig_names.txt" "${outpath}/${SAMPLE}.contigs.fa" > "${outpath2}/${SAMPLE}_RVDB_blastcontigs.fasta"
sed -i 's/--//' "${outpath2}/${SAMPLE}_RVDB_blastcontigs.fasta"
sed -i '/^[[:space:]]*$/d' "${outpath2}/${SAMPLE}_RVDB_blastcontigs.fasta"
sed --posix -i "/^\>/ s/$/"_$SAMPLE"/" "${outpath2}/${SAMPLE}_RVDB_blastcontigs.fasta"
rm "${outpath2}/${SAMPLE}_temp_contig_names.txt"
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/DetectionVirusP1."${SLURM_ARRAY_TASK_ID}".* "${FOLDER}/scripts/sbatch_log"
# ---------------------------------
