#!/bin/bash

# =====================================================
# BactoNano workflow: analysing 16s nanopore sequencing data
# Author: Rodgee Mae Guden

# Copyright © 2025 Rodgee Mae Guden. Marine Biology Research Group, Ghent University
# All rights reserved. This script and its associated workflow are the intellectual property of the author. No part of this code may be reproduced, distributed, or transmitted in any form or by any means, including photocopying, recording, or other electronic or mechanical methods, without the prior written permission of the author, except in the case of brief quotations used in academic or research contexts with proper attribution.

# For permissions, contact: rodgeemae.guden@ugent.be

# =====================================================

# Exit on any error
set -euo pipefail

# Configuration Parameters
THREADS=12
READS_DIR="/home/rguden/Documents/NEMADAPT/ONT_fastq/"
OUTDIR="/home/rguden/Documents/NEMADAPT/output"
KRAKEN_DB="/home/rguden/Documents/databases/silva"
LOG_PATH="$OUTDIR/logs"
LOG_FILE="$OUTDIR/logs/16s_pipeline.log"

# Make the necessary directories
mkdir -p $OUTDIR
mkdir -p $KRAKEN_DB
mkdir -p $LOG_PATH

# To run the script
# bash bactonano_workflow.sh -s qc
# bash bactonano_workflow.sh -s qc,trim

# Function to log messages with timestamps
echo_with_timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Error Handling
trap 'echo_with_timestamp "Pipeline failed at: $(date)"; exit 1' ERR

# ------------------------------------------------------------
# Parse Command-Line Arguments for Specific Steps
# ------------------------------------------------------------
while getopts "s:" opt; do
    case ${opt} in
        s)
            IFS=',' read -r -a STEPS <<< "$OPTARG"  # Split steps by commas
            ;;
        *)
            echo "Usage: $0 -s <step1,step2,...>" | tee -a "$LOG_FILE"
            exit 1
            ;;
    esac
done

if [ ${#STEPS[@]} -eq 0 ]; then
    STEPS=(demultiplex qc trim_adapter trim_chimeric filter postqc classify format krona bracken_setup bracken extract_otu kraken_count bracken_count)
fi

echo_with_timestamp "Steps to run: ${STEPS[*]}"

# ------------------------------------------------------------
# Step Processing Loop
# ------------------------------------------------------------
for STEP in "${STEPS[@]}"; do
    echo_with_timestamp "Starting step: $STEP"
    case "$STEP" in

        concat_fastq)
            echo_with_timestamp "Step: Concatenate FASTQ files per barcode"
            #CONCAT_DIR="$OUTDIR/concatenated/FC3"
            #CONCAT_DIR="$OUTDIR/concatenated/FC2"
            CONCAT_DIR="$OUTDIR/concatenated/FC1"
            mkdir -p "$CONCAT_DIR"

            #BARCODE_PARENT_DIR="$READS_DIR/FC3"
            #BARCODE_PARENT_DIR="$READS_DIR/FC2"
            BARCODE_PARENT_DIR="$READS_DIR/FC1"

            for barcode_dir in "$BARCODE_PARENT_DIR"/barcode*; do
                barcode_name=$(basename "$barcode_dir")
                output_file="$CONCAT_DIR/${barcode_name}.fastq.gz"

                echo_with_timestamp "Concatenating files in $barcode_name..."
                cat "$barcode_dir"/*.fastq.gz > "$output_file"
            done

            echo_with_timestamp "✅ All barcode FASTQ files concatenated to: $CONCAT_DIR"
            ;;

        rename_fastq_FC3)
            echo_with_timestamp "Step: Rename concatenated FASTQ files based on sample names"

            CONCAT_DIR="$OUTDIR/concatenated/FC3"
            RENAME_DIR="$OUTDIR/renamed"
            mkdir -p "$RENAME_DIR"

            # Sample names ordered column-wise
            sample_names=(
                1B23P1 1B23P2 1B23P3 1B23P4 1C20U1 1C20U2 1C20U3 1C20P1
                1C20P2 1C23U1 1C23U2 1C23P1 1C23P2 1D20U1 1D20U2 1D20P1
                1D20P2 1D23U1 1D23U2 1D23U3 1D23P1 1D23P2 blank1 blank2
            )

            for i in "${!sample_names[@]}"; do
                barcode_index=$(printf "%02d" $((i + 1)))
                original_file="$CONCAT_DIR/barcode${barcode_index}.fastq.gz"
                new_file="$RENAME_DIR/${sample_names[$i]}.fastq.gz"

                if [[ -f "$original_file" ]]; then
                    cp "$original_file" "$new_file"
                    echo_with_timestamp "Renamed barcode${barcode_index} → ${sample_names[$i]}"
                else
                    echo_with_timestamp "⚠️ WARNING: $original_file not found!"
                fi
            done

            echo_with_timestamp "✅ Renaming completed. Output in $RENAME_DIR"
            ;;

        rename_fastq_FC1)
            echo_with_timestamp "Step: Rename concatenated FASTQ files"

            CONCAT_DIR="$OUTDIR/concatenated/FC1"
            RENAME_DIR="$OUTDIR/renamed"
            mkdir -p "$RENAME_DIR"

            # Sample names in column-wise order
            sample_names_set3=(
                020U1 020U2 020U3 020U4 020U5 020P1 020P2 020P3
                020P4 020P5 023U1 023U2 023U3 023U4 023U5 023P1
                023P2 023P3 023P4 023P5 1A20U1 1A20U2 1A20U3 1A20U4
            )

            for i in "${!sample_names_set3[@]}"; do
                barcode_index=$(printf "%02d" $((i + 1)))
                original_file="$CONCAT_DIR/barcode${barcode_index}.fastq.gz"
                new_file="$RENAME_DIR/${sample_names_set3[$i]}.fastq.gz"

                if [[ -f "$original_file" ]]; then
                    cp "$original_file" "$new_file"
                    echo_with_timestamp "Renamed barcode${barcode_index} → ${sample_names_set3[$i]}"
                else
                    echo_with_timestamp "⚠️ WARNING: $original_file not found!"
                fi
            done

            echo_with_timestamp "✅ Renaming completed. Output in: $RENAME_DIR"
            ;;

        nanostat)
            echo_with_timestamp "Step: Running NanoStat on renamed FASTQ files"

            INPUT_DIR="$OUTDIR/renamed"
            NANOSTAT_DIR="$OUTDIR/nanostat"
            mkdir -p "$NANOSTAT_DIR"

            COMBINED_TABLE="$NANOSTAT_DIR/combined_nanostat_summary.tsv"
            echo -e "Sample\tNumber of reads\tTotal bases\tMean read length\tMedian read length\tRead length N50\tMean Qscore\tMedian Qscore" > "$COMBINED_TABLE"

            for file in "$INPUT_DIR"/*.fastq.gz; do
                sample_name=$(basename "$file" .fastq.gz)
                output_file="$NANOSTAT_DIR/$sample_name"

                # Run NanoStat
                NanoStat --fastq "$file" --outdir "$NANOSTAT_DIR" --name "$sample_name" &>> "$LOG_FILE"
                echo_with_timestamp "✅ NanoStat completed for $sample_name"

                if [[ ! -f "$output_file" ]]; then
                    echo_with_timestamp "⚠️ WARNING: NanoStat report missing for $sample_name"
                    continue
                fi

                # Extract key values
                num_reads=$(grep -m1 "Number of reads:" "$output_file" | awk -F': ' '{gsub(",", "", $2); print $2}')
                total_bases=$(grep -m1 "Total bases:" "$output_file" | awk -F': ' '{gsub(",", "", $2); print $2}')
                mean_len=$(grep -m1 "Mean read length:" "$output_file" | awk -F': ' '{gsub(",", ".", $2); print $2}')
                median_len=$(grep -m1 "Median read length:" "$output_file" | awk -F': ' '{gsub(",", ".", $2); print $2}')
                n50=$(grep -m1 "Read length N50:" "$output_file" | awk -F': ' '{gsub(",", ".", $2); print $2}')
                mean_qscore=$(grep -m1 "Mean read quality:" "$output_file" | awk -F': ' '{gsub(",", ".", $2); print $2}')
                median_qscore=$(grep -m1 "Median read quality:" "$output_file" | awk -F': ' '{gsub(",", ".", $2); print $2}')

                echo -e "${sample_name}\t${num_reads}\t${total_bases}\t${mean_len}\t${median_len}\t${n50}\t${mean_qscore}\t${median_qscore}" >> "$COMBINED_TABLE"
            done

            echo_with_timestamp "📊 All NanoStat summaries saved to: $NANOSTAT_DIR"
            echo_with_timestamp "📁 Combined summary table saved to: $COMBINED_TABLE"
            ;;

        qc)
            echo_with_timestamp "Step: Quality Control (FastQC & MultiQC)"
            mkdir -p "$OUTDIR/qc"

            INPUT_DIR="$OUTDIR/renamed"
            for file in "$INPUT_DIR"/*.fastq.gz; do
                fastqc -t "$THREADS" -o "$OUTDIR/qc" "$file" &>> "$LOG_FILE"
            done
            multiqc "$OUTDIR/qc" -o "$OUTDIR/qc/multiqc_report" &>> "$LOG_FILE"
            ;;

        trim_adapter)
            echo_with_timestamp "Step: Adapter Removal with Porechop"
            mkdir -p "$OUTDIR/trimmed"

            INPUT_DIR="$OUTDIR/renamed"
            for file in "$INPUT_DIR"/*.fastq.gz; do
                porechop -i "$file" -o "$OUTDIR/trimmed/$(basename "$file" .fastq)_trimmed.fastq.gz" -t "$THREADS" &>> "$LOG_FILE"
            done
            ;;

        trim_chimeric)
            echo_with_timestamp "Step: Detecting and Removing Chimeric Reads"
            mkdir -p "$OUTDIR/trimmed/not_chimeric"

            # Generate overlaps for chimera detection
            for file in "$OUTDIR/trimmed"/*_trimmed.fastq.gz; do
                overlap_paf="$OUTDIR/trimmed/overlap_$(basename "$file").paf"
                minimap2 -x ava-ont -g 500 "$file" "$file" > "$overlap_paf"
            done

            # Detect and remove chimeras
            for file in "$OUTDIR/trimmed"/*_trimmed.fastq.gz; do
                base=$(basename "$file")
                overlap_paf="$OUTDIR/trimmed/overlap_$base.paf"
                report_yacrd="$OUTDIR/trimmed/report_$base.yacrd"
                scrubb_out="$OUTDIR/trimmed/not_chimeric/scrubb_$base.fastq.gz"

                yacrd -i "$overlap_paf" -o "$report_yacrd" -c 4 -n 0.4 scrubb -i "$file" -o "$scrubb_out" &>> "$LOG_FILE"
            done
            ;;

        filter)
            echo_with_timestamp "Step: Filtering with Fastp on non-chimeric reads"
            mkdir -p "$OUTDIR/filtered"
            for file in "$OUTDIR/trimmed/not_chimeric"/scrubb_*.fastq.gz; do
                base=$(basename "$file" .fastq.gz)
                output_file="$OUTDIR/filtered/${base}_filtered.fastq.gz"
                fastp -i "$file" \
                    -o "$output_file" \
                    -q 9 -l 1000 --length_limit 2000 \
                    -h "$OUTDIR/filtered/${base}_fastp.html" \
                    -j "$OUTDIR/filtered/${base}_fastp.json" \
                    -t "$THREADS" &>> "$LOG_FILE"
                echo_with_timestamp "Filtered: $base"
            done
            ;;

        postqc)
            echo_with_timestamp "Step: Post-Filtering QC (FastQC & MultiQC)"
            mkdir -p "$OUTDIR/qc_filtered"
            for file in "$OUTDIR/filtered"/scrubb_*.fastq.gz; do
                fastqc -t "$THREADS" -o "$OUTDIR/qc_filtered" "$file" &>> "$LOG_FILE"
            done
            multiqc "$OUTDIR/qc_filtered" -o "$OUTDIR/qc_filtered/multiqc_report" &>> "$LOG_FILE"
            ;;

        classify)
            echo_with_timestamp "Step: Taxonomic Classification with Kraken2"
            mkdir -p "$OUTDIR/taxonomy"

            # Build and download Kraken database, if needed
            # create db folder (e.g., mkdir -p /home/rguden/Documents/databases/silva)
            # kraken2-build --special silva --db /home/rguden/Documents/databases/silva --threads 4
            # kraken2-build --build --db /home/rguden/Documents/databases/silva --threads 4

            # Run kraken2
            for file in "$OUTDIR/filtered"/scrubb_*.fastq.gz; do
                kraken2 --db "$KRAKEN_DB" \
                        --report "$OUTDIR/taxonomy/$(basename "$file" .fastq)_report.txt" \
                        --output "$OUTDIR/taxonomy/$(basename "$file" .fastq)_kraken.txt" \
                        --threads "$THREADS" "$file" &>> "$LOG_FILE"
            done
            ;;

        krona)
            echo_with_timestamp "Step: Visualization with Krona (Multiple Samples)"
            mkdir -p "$OUTDIR/visualization"

            # Load krona taxonomy database
            # ktUpdateTaxonomy.sh

            # Check if there are any formatted Krona input files
            if ls "$OUTDIR/taxonomy"/*_report.txt 1> /dev/null 2>&1; then
                # Generate a single Krona chart combining multiple samples
                echo_with_timestamp "Generating combined Krona chart from multiple samples"

                # Run ktImportTaxonomy with multiple input files
                ktImportTaxonomy -t 5 -m 3 -o "$OUTDIR/visualization/multi-krona.html" "$OUTDIR/taxonomy"/*_report.txt &>> "$LOG_FILE"

                # Check if the multi-sample Krona file was successfully created
                if [[ -f "$OUTDIR/visualization/multi-krona.html" ]]; then
                    echo_with_timestamp "Generated combined Krona chart: $OUTDIR/visualization/multi-krona.html"
                else
                    echo_with_timestamp "Error: Failed to generate combined Krona chart"
                fi
            else
                echo_with_timestamp "Warning: No Krona input files found in $OUTDIR/krona"
            fi
            ;;

        bracken_setup)
            echo_with_timestamp "Step: Setting up Bracken database"

            # Generate the Bracken database file
            echo_with_timestamp "Generating Bracken database file"
            bracken-build -d "$KRAKEN_DB" -t "$THREADS" -k 35 -l 1500 &>> "$LOG_FILE"

            echo_with_timestamp "Bracken database setup completed."
            ;;

        bracken)
            echo_with_timestamp "Step: Estimating species abundance with Bracken"
            mkdir -p "$OUTDIR/bracken"
            for file in "$OUTDIR/taxonomy"/*_report.txt; do
                bracken -d "$KRAKEN_DB" -i "$file" -o "$OUTDIR/bracken/$(basename "$file" _report.txt)_bracken_genus.txt" -r 1500 -l G &>> "$LOG_FILE"
                #bracken -d "$KRAKEN_DB" -i "$file" -o "$OUTDIR/bracken/$(basename "$file" _report.txt)_bracken_genus.txt" -r 1500 -l S &>> "$LOG_FILE"
                echo_with_timestamp "Bracken abundance estimation completed for $file"
            done
            ;;

         extract_otu)
            echo_with_timestamp "Step: Extract OTU Tables from Bracken Reports (All Levels)"
            mkdir -p "$OUTDIR/otu_tables"

            # Define ranks and corresponding filename stems
            declare -A RANK_NAMES=( ["S"]="species" ["G"]="genus" ["F"]="family" ["O"]="order" ["C"]="class" ["P"]="phylum" )
            for LEVEL in S G F O C P; do
                RANK_NAME=${RANK_NAMES[$LEVEL]}
                OTU_TABLE="$OUTDIR/otu_tables/otu_table_${RANK_NAME}_bracken.tsv"
                TMP="$OUTDIR/otu_tables/tmp_otu_bracken_${RANK_NAME}.tsv"

                > "$TMP"

                # Process all Bracken files for this rank
                for file in "$OUTDIR/bracken"/*_bracken_${RANK_NAME}.txt; do
                    if [[ -f "$file" ]]; then
                        sample_name=$(basename "$file" .txt)
                        echo_with_timestamp "Processing Bracken OTU ($RANK_NAME) from file: $file"

                        tail -n +2 "$file" | awk -F'\t' -v level="$LEVEL" -v sname="$sample_name" '
                            $3 == level {
                                # Clean taxon name
                                gsub(/[[:space:]]+/, "_", $1)
                                gsub(/_+/, "_", $1)
                                gsub(/^_+|_+$/, "", $1)
                                print $1 "\t" $3 "\t" $6 "\t" sname
                            }
                        ' >> "$TMP"
                    fi
                done

                # Skip if no rows in TMP (i.e., no files for this rank)
                if [[ ! -s "$TMP" ]]; then
                    echo_with_timestamp "No Bracken files found for $RANK_NAME. Skipping OTU table."
                    continue
                fi

                # Collect unique taxon names and sample names
                awk -F'\t' '{print $1}' "$TMP" | sort -u > "$OUTDIR/otu_tables/unique_${RANK_NAME}_names.txt"
                awk -F'\t' '{print $4}' "$TMP" | sort -u > "$OUTDIR/otu_tables/unique_${RANK_NAME}_samples.txt"

                # Build OTU table header (TSV)
                {
                    printf "otu"
                    while read sample; do printf "\t%s" "$sample"; done < "$OUTDIR/otu_tables/unique_${RANK_NAME}_samples.txt"
                    echo

                    while read taxon; do
                        printf "%s" "$taxon"
                        while read sample; do
                            value=$(awk -F'\t' -v g="$taxon" -v s="$sample" '$1 == g && $4 == s {print $3}' "$TMP")
                            printf "\t%s" "${value:-0}"
                        done < "$OUTDIR/otu_tables/unique_${RANK_NAME}_samples.txt"
                        echo
                    done < "$OUTDIR/otu_tables/unique_${RANK_NAME}_names.txt"
                } > "$OTU_TABLE"

                # Clean up temp files for this rank
                rm "$TMP" "$OUTDIR/otu_tables/unique_${RANK_NAME}_names.txt" "$OUTDIR/otu_tables/unique_${RANK_NAME}_samples.txt"

                echo_with_timestamp "OTU table generated for $RANK_NAME: $OTU_TABLE"
            done

            echo_with_timestamp "All OTU tables generated successfully for Bracken"
            ;;

        *)
            echo_with_timestamp "Unknown step: $STEP" &>> "$LOG_FILE"
            ;;
    esac
    echo_with_timestamp "Completed step: $STEP"
done
