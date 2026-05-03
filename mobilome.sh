#!/bin/bash

# === Configuration ===
INPUT_DIR="Genomes"
OUT_DIR="mobilome_results"
THREADS=32
#CHECKV_DB="/home/checkv-db-v1.5"  # Modify to your CheckV DB location

# === Create Output Directory ===
mkdir -p "$OUT_DIR"

# === Initialize Summary Files ===
echo -e "Genome\tPlasmid_ID\tReplicon\tRelaxase" > "$OUT_DIR/mob_suite_summary.tsv"
#echo -e "Genome\tCheckV_Contig\tCheckV_Category\tGenome_Size\tCompleteness" > "$OUT_DIR/checkv_summary.tsv"
echo -e "Genome\tIntegron\tType\tGene_Cassette_Info" > "$OUT_DIR/integron_summary.tsv"
echo -e "Genome\tIS_Element\tFamily\tStart\tEnd" > "$OUT_DIR/isescan_summary.tsv"

# === Loop through all genome FASTA files ===
for genome in "$INPUT_DIR"/*.fasta; do
    base=$(basename "$genome" .fasta)
    echo "🔬 Processing $base"

    #### MOB-suite Analysis ####
    mob_recon --infile "$genome" --outdir "$OUT_DIR/${base}_mob" --num_threads "$THREADS"
    mobtyper_file="$OUT_DIR/${base}_mob/mobtyper_results.txt"

    if [ -f "$mobtyper_file" ]; then
        awk -F '\t' -v gname="$base" 'NR>1 {
            split($1, s, ":");  # s[2] is the plasmid ID
            print gname"\t"s[2]"\t"$6"\t"$8
        }' "$mobtyper_file" >> "$OUT_DIR/mob_suite_summary.tsv"
    else
        echo "⚠️  MOB-suite results not found for $base, skipping MOB-suite summary entry."
    fi

    #### CheckV Analysis ####
    #mkdir -p "$OUT_DIR/${base}_checkv"
    #checkv end_to_end "$genome" "$OUT_DIR/${base}_checkv" -t "$THREADS" -d "$CHECKV_DB"
    #if [ -f "$OUT_DIR/${base}_checkv/quality_summary.tsv" ]; then
    #    awk -F '\t' -v gname="$base" 'NR>1 {
    #        print gname"\t"$1"\t"$2"\t"$3"\t"$4
    #     }' "$OUT_DIR/${base}_checkv/quality_summary.tsv" >> "$OUT_DIR/checkv_summary.tsv"
    #fi


    #### IntegronFinder Analysis ####
    integron_finder "$genome" --cpu "$THREADS" --func-annot --local-max --outdir "$OUT_DIR/${base}_integrons"
    integron_file="$OUT_DIR/${base}_integrons/Results_Integron_Finder_${base}/${base}.integrons"
    if [ -f "$integron_file" ]; then
        awk -F '\t' -v gname="$base" 'NR>1 {
            print gname"\t"$1"\t"$2"\t"$4"\t"$5"\t"$8"\t"$9
        }' "$integron_file" >> "$OUT_DIR/integron_summary.tsv"
    fi

    #### ISEScan Analysis ####
    isescan.py --seqfile "$genome" --output "$OUT_DIR/${base}_isescan" --nthread "$THREADS"
    isescan_file="$OUT_DIR/${base}_isescan/Genomes/${base}.fasta.tsv"
    if [ -f "$isescan_file" ]; then
        awk -F '\t' -v gname="$base" 'NR>1 {
            print gname"\t"$1"\t"$2"\t"$4"\t"$5"\t"$18"\t"$20"\t"$22
        }' "$isescan_file" >> "$OUT_DIR/isescan_summary.tsv"
    fi
    
    echo "✅ Finished $base"
done
echo "🎯 Mobilome analysis completed for all samples."

