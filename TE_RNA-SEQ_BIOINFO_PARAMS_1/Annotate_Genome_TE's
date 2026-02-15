#!/bin/bash
#SBATCH --job-name=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=rm%j.out
#SBATCH --error=rm%j.err

# ============================================================================
# TE Annotation Pipeline Using Custom Consensus Library
# ============================================================================

#genome.fa is the flybase DM6 (6.64) and te_consensus.fa is the Bergman labs D.mel TE concensus sequences

set -e

GENOME="genome.fa"
TE_LIBRARY="te_consensus.fa"
THREADS=8
OUTPUT_DIR="TE_repeatmasker"

mkdir -p ${OUTPUT_DIR}

echo "============================================================================"
echo "Running RepeatMasker with custom TE consensus library"
echo "============================================================================"
echo ""
echo "Genome: ${GENOME}"
echo "TE library: ${TE_LIBRARY}"
echo "Threads: ${THREADS}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check files exist
if [ ! -f "${GENOME}" ]; then
    echo "ERROR: Genome file ${GENOME} not found!"
    exit 1
fi

if [ ! -f "${TE_LIBRARY}" ]; then
    echo "ERROR: TE library ${TE_LIBRARY} not found!"
    exit 1
fi

if [ ! -f "${GENOME}.fai" ]; then
    echo "Creating genome index..."
    samtools faidx ${GENOME}
fi

# ============================================================================
# Step 1: Run RepeatMasker
# ============================================================================

echo "Step 1: Running RepeatMasker..."
echo "============================================================================"

RepeatMasker \
    -lib ${TE_LIBRARY} \
    -pa ${THREADS} \
    -s \
    -gff \
    -nolow \
    -dir ${OUTPUT_DIR} \
    ${GENOME}

echo ""
echo "RepeatMasker complete!"
echo ""

# ============================================================================
# Step 2: Convert to GTF
# ============================================================================

echo "Step 2: Converting RepeatMasker output to GTF..."
echo "============================================================================"

OUTPUT_FILE="${OUTPUT_DIR}/$(basename ${GENOME})"

python3 repeatmasker_to_gtf.py \
    ${OUTPUT_FILE}.out \
    ${OUTPUT_DIR}/TE_annotations_complete.gtf

echo ""
echo "GTF conversion complete!"
echo ""

# ============================================================================
# Step 3: Generate Statistics
# ============================================================================

echo "Step 3: Generating statistics..."
echo "============================================================================"

GTF_FILE="${OUTPUT_DIR}/TE_annotations_complete.gtf"

# Count total TEs
TOTAL_TES=$(grep -v "^#" ${GTF_FILE} | wc -l)
echo "Total TE insertions: $(printf "%'d" ${TOTAL_TES})"

# Count unique families
FAMILIES=$(grep -v "^#" ${GTF_FILE} | \
    awk -F'family_id "' '{print $2}' | \
    awk -F'"' '{print $1}' | \
    sort -u | wc -l)
echo "Unique TE families: ${FAMILIES}"

# Calculate genome coverage
TOTAL_BP=$(grep -v "^#" ${GTF_FILE} | \
    awk '{sum += $5 - $4 + 1} END {print sum}')
GENOME_SIZE=$(awk '{sum += $2} END {print sum}' ${GENOME}.fai)
COVERAGE=$(awk -v tb=${TOTAL_BP} -v gs=${GENOME_SIZE} 'BEGIN {printf "%.2f", (tb/gs)*100}')

echo "TE coverage: ${COVERAGE}%"
echo "Total TE bp: $(printf "%'d" ${TOTAL_BP})"
echo "Genome size: $(printf "%'d" ${GENOME_SIZE})"
echo ""

# Count by chromosome
echo "TE count by chromosome:"
echo "-----------------------"
grep -v "^#" ${GTF_FILE} | \
    awk '{print $1}' | \
    sort | uniq -c | sort -k2 | \
    awk '{printf "  %s: %'"'"'d\n", $2, $1}'

echo ""
echo "Top 20 TE families:"
echo "-------------------"
grep -v "^#" ${GTF_FILE} | \
    awk -F'family_id "' '{print $2}' | \
    awk -F'"' '{print $1}' | \
    sort | uniq -c | sort -rn | head -20 | \
    awk '{printf "  %-20s %'"'"'d\n", $2, $1}'

# ============================================================================
# Step 4: Compare with original FlyBase annotation
# ============================================================================

if [ -f "TE_annotations.gtf" ]; then
    echo ""
    echo "Step 4: Comparing with original FlyBase annotation..."
    echo "============================================================================"
    
    OLD_COUNT=$(grep -v "^#" TE_annotations.gtf | wc -l)
    OLD_BP=$(grep -v "^#" TE_annotations.gtf | awk '{sum += $5 - $4 + 1} END {print sum}')
    OLD_COVERAGE=$(awk -v tb=${OLD_BP} -v gs=${GENOME_SIZE} 'BEGIN {printf "%.2f", (tb/gs)*100}')
    
    echo "Original FlyBase GTF:"
    echo "  Insertions: $(printf "%'d" ${OLD_COUNT})"
    echo "  Coverage: ${OLD_COVERAGE}%"
    echo ""
    echo "New RepeatMasker GTF:"
    echo "  Insertions: $(printf "%'d" ${TOTAL_TES})"
    echo "  Coverage: ${COVERAGE}%"
    echo ""
    
    INCREASE=$(awk -v new=${TOTAL_TES} -v old=${OLD_COUNT} 'BEGIN {printf "%.1f", ((new-old)/old)*100}')
    echo "Improvement: +${INCREASE}% more TE insertions detected"
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "============================================================================"
echo "Pipeline Complete!"
echo "============================================================================"
echo ""
echo "Output files:"
echo "  RepeatMasker output: ${OUTPUT_FILE}.out"
echo "  GTF annotation: ${GTF_FILE}"
echo "  Masked genome: ${OUTPUT_FILE}.masked"
echo "  Summary table: ${OUTPUT_FILE}.tbl"
echo ""
echo "To visualize with circos plot:"
echo "  # Update the script to use the new GTF"
echo "  sed -i 's|GTF_FILE <- \"TE_annotations.gtf\"|GTF_FILE <- \"${GTF_FILE}\"|' te_circos_fixed.R"
echo "  Rscript te_circos_fixed.R"
echo ""
echo "============================================================================"
