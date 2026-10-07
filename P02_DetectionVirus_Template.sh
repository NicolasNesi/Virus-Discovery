#!/bin/bash
#SBATCH -J DetectionVirusP2                     # Job name 
#SBATCH -o DetectionVirusP2.%a.%A.out           # File to which stdout will be written
#SBATCH -e DetectionVirusP2.%a.%A.err           # File to which stderr will be written
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
# diamond=2.1.8
# seqkit=2.6.1
# ---------------------------------

# ---------------------------------
# Variables paths
# ---------------------------------
GROUP="${1}"
PROJECT="${2}"
FOLDER="/home/2019013/Data/$GROUP/$PROJECT"
inpath="${FOLDER}/contigs/final_contigs"
outpath="${FOLDER}/blast_results"
tempdir="${FOLDER}/blast_results"
rawreads="${FOLDER}/raw_reads"
readcount="${FOLDER}/read_count"
abundance="${FOLDER}/abundance/final_abundance"

# Database for BLAST
dbnt="/home/2019013/Databases/BLAST_NCBI/nt"
dbnr="home/2019013/Databases/BLAST_NCBI/nr"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
# Concatenate all rdrp and RVDB Blast Contigs
# ---------------------------------
cat "${outpath2}"/*_rdrp_blastcontigs.fasta > "${inpath}/Combined_rdrp_blastcontigs_forNT.fasta"
cat "${outpath2}"/*_RVDB_blastcontigs.fasta > "${inpath}/Combined_RVDB_blastcontigs_forNT.fasta"
cat "${inpath}/Combined_rdrp_blastcontigs_forNT.fasta" "${inpath}/Combined_RVDB_blastcontigs_forNT.fasta" > "${inpath}/Combined-rdrp-RVDB_blastcontigs_forNT.fasta"
# ---------------------------------

# ---------------------------------
# BLAST blastn versus nt databases
# ---------------------------------
echo "Starting Blastn of rdrp-RVDB combined contigs against the database nt"

blastn \
-query "${inpath}/Combined-rdrp-RVDB_blastcontigs_forNT.fasta" \
-db "${dbnt}" \
-out "${outpath}/Combined-rdrp-RVDB_nt_blastn_results.txt" \
-max_target_seqs 10 \
-num_threads "${THREADS}" \
-mt_mode 1 \
-evalue 1E-10 \
-subject_besthit \
-outfmt '6 qseqid qlen sacc salltitles staxids pident length evalue'
# ---------------------------------

# ---------------------------------
# Extract contigs from BLAST to fasta
# ---------------------------------
grep -i ".*" "${outpath}/Combined-rdrp-RVDB_nt_blastn_results.txt" | cut -f1 | sort | uniq > "${outpath}/Combined-rdrp-RVDB_temp_nt_contig_names.txt" 
grep -A1 -I -Ff "${outpath}/Combined-rdrp-RVDB_temp_nt_contig_names.txt" "${inpath}/Combined-rdrp-RVDB_blastcontigs_forNT.fasta" > "${outpath}/Combined-rdrp-RVDB_nt_blastcontigs.fasta"
sed -i 's/--//' "${outpath}/Combined-rdrp-RVDB_nt_blastcontigs.fasta"
sed -i '/^[[:space:]]*$/d' "${outpath}/Combined-rdrp-RVDB_nt_blastcontigs.fasta"
sed --posix -i "/^\>/ s/$/"_$SAMPLE"/" "${outpath}/Combined-rdrp-RVDB_nt_blastcontigs.fasta"
rm "${outpath}/Combined-rdrp-RVDB_temp_nt_contig_names.txt"
# ---------------------------------

# ---------------------------------
# Concatenate all rdrp and RVDB Blast Contigs
# ---------------------------------
cat "${outpath2}"/*_rdrp_blastcontigs.fasta > "${inpath}/Combined_rdrp_blastcontigs_forNR.fasta"
cat "${outpath2}"/*_RVDB_blastcontigs.fasta > "${inpath}/Combined_RVDB_blastcontigs_forNR.fasta"
cat "${inpath}/Combined_rdrp_blastcontigs_forNR.fasta" "${inpath}/Combined_RVDB_blastcontigs_forNR.fasta" > "${inpath}/Combined-rdrp-RVDB_blastcontigs_forNR.fasta"
# ---------------------------------

# ---------------------------------
# Diamond blastx versus nr databases
# ---------------------------------
echo "Starting Blastx of rdrp-RVDB combined contigs against the database nr"

diamond blastx \
--query "${inpath}/Combined-rdrp-RVDB_blastcontigs_forNR.fasta" \
--db "${dbnr}" \
--tmpdir "${tempdir}" \
--out "${outpath}/Combined-rdrp-RVDB_nr_blastx_results.txt" \
--evalue 1E-4 \
-c1 \
--threads "${THREADS}" \
--max-target-seqs 10 \
--block-size 14 \
--outfmt 6 qseqid qlen sseqid stitle staxids pident length evalue \
--more-sensitive
# ---------------------------------

# ---------------------------------
# Extract contigs from BLAST to fasta
# ---------------------------------
grep -i ".*" "${outpath}/Combined-rdrp-RVDB_nr_blastx_results.txt" | cut -f1 | sort | uniq > "${outpath}/Combined-rdrp-RVDB_temp_nr_contig_names.txt" 
grep -A1 -I -Ff "${outpath}/Combined-rdrp-RVDB_temp_nr_contig_names.txt" "${inpath}/Combined-rdrp-RVDB_blastcontigs_forNR.fasta" > "${outpath}/Combined-rdrp-RVDB_nr_blastcontigs.fasta"
sed -i 's/--//' "${outpath}/Combined-rdrp-RVDB_nr_blastcontigs.fasta"
sed -i '/^[[:space:]]*$/d' "${outpath}/Combined-rdrp-RVDB_nr_blastcontigs.fasta"
sed --posix -i "/^\>/ s/$/"_$SAMPLE"/" "${outpath}/Combined-rdrp-RVDB_nr_blastcontigs.fasta"
rm "${outpath}/Combined-rdrp-RVDB_temp_nr_contig_names.txt"
# ---------------------------------

# ---------------------------------
# Count number of reads
# ---------------------------------
cd "${rawreads}/"
seqkit stats --threads 96 -Ta *.fastq.gz | csvtk cut -t -f 1,4 | csvtk del-header > "${reascount}/${PROJECT}_accessions_reads"
# ---------------------------------

# ---------------------------------
# Create summary talbe of Blast outputs
# ---------------------------------
mkdir -p "${FOLDER}/blast_results/summary_table_creation"

abundance_files=("${abundance}"/*_RSEM.isoforms.results)
rdrp_blast_files=("${FOLDER}"/blast_results/*_rdrp_blastx_results.txt)
rvdb_blast_files=("${FOLDER}"/blast_results/*_RVDB_blastx_results.txt)
rdrp_blastcontigs_file=("${FOLDER}"/blast_results/*_rdrp_blastcontigs.fasta)
rvdb_blastcontigs_file=("${FOLDER}"/blast_results/*_RVDB_blastcontigs.fasta)

echo "Number of Abundance files loaded: ${#abundance_files[@]}"

for i in "${abundance_files[@]}"; do
   id=$(basename "$i" | awk '{split($0,a,"_RSEM.isoforms.results"); print a[1]}')
   awk -F'\t' -vOFS='\t' '{ $1 = $1 "_" id }1' id="$id" "$i" > "${abundance}"/"${id}"_anno_RSEM.isoforms.results
done

cat "${abundance}"/*_anno_RSEM.isoforms.results > "${FOLDER}"/blast_results/summary_table_creation/combined_abundance_table.txt
rm "${abundance}"/*_anno_RSEM.isoforms.results 

echo "Number of RdRp Blast files loaded: ${#rdrp_blast_files[@]}"

for i in "${rdrp_blast_files[@]}"; do
    id=$(basename "$i" | awk '{split($0,a,"_rdrp_blastx_results.txt"); print a[1]}')                                  
    awk -F'\t' -vOFS='\t' '{ $1 = $1 "_" id }1' id="$id" "$i" > "${FOLDER}"/blast_results/"$id"_anno_rdrp_blastx_results.txt
done

cat "${FOLDER}"/blast_results/*_anno_rdrp_blastx_results.txt > "${FOLDER}"/blast_results/summary_table_creation/combined_rdrp_blastx_results.txt 
rm "${FOLDER}"/blast_results/*_anno_rdrp_blastx_results.txt 

echo "Number of RVDB Blast files loaded: ${#rvdb_blast_files[@]}"

for i in "${rvdb_blast_files[@]}"; do
   id=$(basename "$i" | awk '{split($0,a,"_RVDB_blastx_results.txt"); print a[1]}')
   awk -F'\t' -vOFS='\t' '{ $1 = $1 "_" id }1' id="$id" "$i" > "${FOLDER}"/blast_results/"$id"_anno_RVDB_blastx_results.txt
done

cat "${FOLDER}"/blast_results/*_anno_RVDB_blastx_results.txt > "${FOLDER}"/blast_results/summary_table_creation/combined_RVDB_blastx_results.txt
rm "${FOLDER}"/blast_results/*_anno_RVDB_blastx_results.txt

echo "Creating a combined virus contigs file"

cat "${rdrp_blastcontigs_file[@]}" > "${FOLDER}"/blast_results/summary_table_creation/combined_rdrp_blastcontigs.fasta
cat "${rvdb_blastcontigs_file[@]}" > "${FOLDER}"/blast_results/summary_table_creation/combined_RVDB_blastcontigs.fasta
cat "${FOLDER}"/blast_results/summary_table_creation/combined_rdrp_blastcontigs.fasta > "${FOLDER}"/blast_results/summary_table_creation/combined_contigs.fasta
cat "${FOLDER}"/blast_results/summary_table_creation/combined_RVDB_blastcontigs.fasta >> "${FOLDER}"/blast_results/summary_table_creation/combined_contigs.fasta
rm "${FOLDER}"/blast_results/summary_table_creation/combined_rdrp_blastcontigs.fasta "${FOLDER}"/blast_results/summary_table_creation/combined_RVDB_blastcontigs.fasta


file_of_accessions_without_path=$(basename "$file_of_accessions")
file_of_accessions_name=$(echo "${file_of_accessions_without_path%.*}")
nr_results="${FOLDER}/blast_results/${SAMPLE}_nr_blastx_results.txt"
nt_results="${FOLDER}/blast_results/${SAMPLE}_nt_blastn_results.txt"
nr_blastcontigs="${FOLDER}/blast_results/${SAMPLE}_nr_blastcontigs.fasta"
nt_blastcontigs="${FOLDER}/blast_results/${SAMPLE}_nt_blastcontigs.fasta"

echo "Process RVDB taxa information"

cut -f3 "${FOLDER}"/blast_results/summary_table_creation/combined_RVDB_blastx_results.txt | cut -f3 -d"|" | sort -u \
> "${FOLDER}"/blast_results/summary_table_creation//temp_joint_RVDB_blast_table_taxids

grep -Ff "${FOLDER}"/blast_results/summary_table_creation//temp_joint_RVDB_blast_table_taxids "$rvdb_accession2taxid" \
> "${FOLDER}"/blast_results/summary_table_creation/temp_joint_RVDB_blast_table_taxids_filtered

cut -f4 "${FOLDER}"/blast_results/summary_table_creation/combined_RVDB_blastx_results.txt | sort -u \
> "${FOLDER}"/blast_results/summary_table_creation/temp_joint_RVDB_blast_table_headers

awk -F '|' '{source=$2 protein_accession=$3 source2=$4 nucl_accession=$5 genomic_region=$6 
            organism=gensub(/.*\[|\].*/, "", "g", genomic_region)
			gsub(/\[.*]/, "", genomic_region)
			gsub(/^[ \t]+[ \t]+$/, "",genomic_region)
			original_line=$0
			gsub(/\|/, "%", original_line)
			print original_line "|" source "|" protein_accession "|" source2 "|" nucl_accession "|" genomic_region "|" organism
			}' "${FOLDER}"/blast_results/summary_table_creation/temp_joint_RVDB_blast_table_headers > "${FOLDER}"/blast_results/summary_table_creation/temp_joint_RVDB_blast_table_headers_parsed
			
echo "Join the two tables"

   join_tables() {
        # Check if both input files exist
        if [ ! -f "$1" ] || [ ! -f "$2" ]; then
            echo "Input files not found."
            return 1
        fi

        # Sort the first input file based on the third column
        sorted_table1=$(mktemp)
        sort -t '|' -k 3,3 "$1" >"$sorted_table1"

        # Sort the second input file based on the second column
        sorted_table2=$(mktemp)
        sort -t $'\t' -k 2,2 "$2" >"$sorted_table2"
        sed -i 's/\t/\|/g' "$sorted_table2"

        # Join the sorted tables using the third column for table 1 and the second column for table 2
        join -t '|' -1 3 -2 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.3 "$sorted_table1" "$sorted_table2"

        # Clean up temporary files
        rm "$sorted_table1" "$sorted_table2"
    }

join_tables "${FOLDER}"/blast_results/summary_table_creation/temp_joint_RVDB_blast_table_headers_parsed "${FOLDER}"/blast_results/summary_table_creation/temp_joint_RVDB_blast_table_taxids_filtered \
> "${FOLDER}"/blast_results/summary_table_creation/temp_joint_RVDB_taxa_blast_table.txt

Rscript "${FOLDER}"/Scripts/create_blast_joint_table.R \
--nr "${FOLDER}/blast_results/summary_table_creation/$nr_results_file" \
--nt "${FOLDER}/blast_results/summary_table_creation/$nt_results_file" \
--rdrp "${FOLDER}/blast_results/summary_table_creation/combined_rdrp_blastx_results.txt" \
--rvdb "${FOLDER}/blast_results/summary_table_creation/combined_RVDB_blastx_results.txt" \
--abundance "${FOLDER}/blast_results/summary_table_creation/combined_abundance_table.txt" \
--readcounts "${reascount}/${PROJECT}_accessions_reads" \
--output "${FOLDER}/blast_results/summary_table_creation/temp_joint_blast_table" \
--rdrp_tax "$taxonomy"/RdRp-scan_0.90.info \
--rvdb_tax "${FOLDER}/blast_results/summary_table_creation/temp_joint_RVDB_taxa_blast_table.txt" \
--multi_lib

echo "Getting lineage for taxids"

current_month="$(date +%b-%Y)"

taxonkit \
lineage -c "${FOLDER}"/blast_results/summary_table_creation/temp_joint_blast_table_taxids \
--data-dir /scratch/VELAB/Databases/Blast/taxdmp."${current_month}"/ | awk '$2>0' | cut -f 2- | \
taxonkit reformat \
--output-ambiguous-result \
--data-dir /scratch/VELAB/Databases/Blast/taxdmp."${current_month}"/ \
-I 1 \
-r "Unassigned" \
-R "missing_taxid" \
--fill-miss-rank -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" | csvtk add-header -t -n "taxid,lineage,kingdom,phylum,class,order,family,genus,species" \
> "${FOLDER}/blast_results/summary_table_creation/temp_lineage_table"


# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/DetectionVirusP2.* "${FOLDER}/scripts/sbatch_log"
# ---------------------------------





