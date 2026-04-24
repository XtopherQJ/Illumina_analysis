#!/bin/bash

# Enable Conda environments

eval "$(conda shell.bash hook)"

echo "Activating the alignment environment..."
conda activate alingment
# === Step 0: Setup ===

RAW_DIR="data"
TRIMMED_DIR="trimmed_reads"
ASSEMBLY_DIR="spades_output"
DRAFT_DIR="draft_assemblies"
PILON_DIR="pilon_output"
QUAST_DIR="quast_results"
ANNOT_DIR="annotation"
ADAPTERS="TruSeq3-PE.fa"
THREADS=8
REFERENCE="reference.fasta"  # Optional for QUAST

mkdir -p "$TRIMMED_DIR" "$ASSEMBLY_DIR" "$DRAFT_DIR" "$PILON_DIR" "$QUAST_DIR" "$ANNOT_DIR"

# === Step 1: Quality Control with FastQC and MultiQC ===
echo "🔍 Running FastQC..."
fastqc "$RAW_DIR"/*.fastq.gz -o fastqc_reports -t $THREADS
multiqc fastqc_reports -o fastqc_reports

# === Step 2: Trimming reads with Trimmomatic ===
echo "✂️ Trimming reads..."
for sample in $(ls "$RAW_DIR"/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | xargs -n1 basename); do
    R1="${RAW_DIR}/${sample}_R1.fastq.gz"
    R2="${RAW_DIR}/${sample}_R2.fastq.gz"
    OUT1="${TRIMMED_DIR}/${sample}_R1_paired.fastq.gz"
    OUT2="${TRIMMED_DIR}/${sample}_R2_paired.fastq.gz"
    UNP1="${TRIMMED_DIR}/${sample}_R1_unpaired.fastq.gz"
    UNP2="${TRIMMED_DIR}/${sample}_R2_unpaired.fastq.gz"

    trimmomatic PE -threads $THREADS -phred33 "$R1" "$R2" "$OUT1" "$UNP1" "$OUT2" "$UNP2" ILLUMINACLIP:$ADAPTERS:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
done

# === Step 3: Genome assembly with SPAdes ===
echo "Activating the SPAdes environment..."
conda activate spades
echo "🧬 Running SPAdes assembly..."
for sample in $(ls "$TRIMMED_DIR"/*_R1_paired.fastq.gz | sed 's/_R1_paired.fastq.gz//' | xargs -n1 basename); do
    spades.py -1 "${TRIMMED_DIR}/${sample}_R1_paired.fastq.gz" \
              -2 "${TRIMMED_DIR}/${sample}_R2_paired.fastq.gz" \
              -o "${ASSEMBLY_DIR}/${sample}_spades" -t $THREADS --careful
    cp "${ASSEMBLY_DIR}/${sample}_spades/contigs.fasta" "${DRAFT_DIR}/${sample}.fasta"
done

# === Step 4: Polishing with Pilon ===
echo "Activating the alignment environment..."
conda activate alingment
echo "🔧 Polishing assemblies with Pilon..."
for draft in "$DRAFT_DIR"/*.fasta; do
    sample=$(basename "$draft" .fasta)
    R1="${TRIMMED_DIR}/${sample}_R1_paired.fastq.gz"
    R2="${TRIMMED_DIR}/${sample}_R2_paired.fastq.gz"
    bwa index "$draft"
    bwa mem -t $THREADS "$draft" "$R1" "$R2" | samtools view -Sb - > "${PILON_DIR}/${sample}.bam"
    samtools sort -o "${PILON_DIR}/${sample}_sorted.bam" "${PILON_DIR}/${sample}.bam"
    samtools index "${PILON_DIR}/${sample}_sorted.bam"
    pilon --genome "$draft" --frags "${PILON_DIR}/${sample}_sorted.bam" --output "${sample}_polished" --outdir "$PILON_DIR" --threads $THREADS
done

echo "Activating the QUAST environment..."
conda activate quast

# === Step 5: Quality Assessment with QUAST ===
echo "📊 Running QUAST..."
quast.py "$DRAFT_DIR"/*.fasta "$PILON_DIR"/*_polished.fasta -o "$QUAST_DIR" -t $THREADS -r "$REFERENCE" || true

echo "Activating the alignment environment..."
conda activate alingment

# === Step 6: Annotation with Prokka ===
echo "🧬 Annotating genomes with Prokka..."
for genome in "$PILON_DIR"/*.fasta; do
    sample=$(basename "$genome" .fasta)
    prokka --cpus $THREADS --kingdom Bacteria --locustag "$sample" --addgenes --centre X --compliant --force --outdir "${ANNOT_DIR}/${sample}" "$genome"
done

echo "✅ Pipeline completed successfully."
