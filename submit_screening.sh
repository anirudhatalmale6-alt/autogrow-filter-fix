#!/bin/bash
#=====================================================================
# Virtual Screening Submission Script for SDSC Expanse
# Screens NP Atlas and/or COCONUT databases using AutoDock Vina
# Runs parallel docking processes on a single node
#
# Usage:
#   bash submit_screening.sh npatlas    # Screen NP Atlas only
#   bash submit_screening.sh coconut    # Screen COCONUT only
#   bash submit_screening.sh both       # Screen both databases
#=====================================================================

# ── Configuration ──
PROJECT_DIR=~/final_project
RECEPTOR="${PROJECT_DIR}/receptors/2ZSH.pdbqt"
DB_DIR="${PROJECT_DIR}/np_databases"
RESULTS_DIR="${PROJECT_DIR}/screening_results"
POSES_DIR="${RESULTS_DIR}/top_poses"
SCREEN_SCRIPT="${PROJECT_DIR}/screen_database.py"

# Docking parameters
CENTER_X=55.0
CENTER_Y=59.5
CENTER_Z=35.5
SIZE_X=21.0
SIZE_Y=21.0
SIZE_Z=22.5
EXHAUSTIVENESS=20

# Slurm settings (matching your normal AutoGrow4 setup)
PARTITION="shared"
ACCOUNT="--account=iit135"
NTASKS=100
MEM="200G"
TIME="48:00:00"

# Parallel workers within the node
PARALLEL_WORKERS=50
# Compounds per worker chunk
CHUNK_SIZE=200

mkdir -p "${RESULTS_DIR}" "${POSES_DIR}" "${RESULTS_DIR}/logs"

# ── Functions ──

count_lines() {
    wc -l < "$1" | tr -d ' '
}

submit_database() {
    local DB_NAME="$1"
    local SMI_FILE="${DB_DIR}/${DB_NAME}.smi"

    if [ ! -f "$SMI_FILE" ]; then
        echo "ERROR: ${SMI_FILE} not found. Run prepare_np_databases.py first."
        return 1
    fi

    local TOTAL=$(count_lines "$SMI_FILE")
    local NUM_CHUNKS=$(( (TOTAL + CHUNK_SIZE - 1) / CHUNK_SIZE ))

    echo "Database: ${DB_NAME}"
    echo "  Compounds: ${TOTAL}"
    echo "  Chunk size: ${CHUNK_SIZE}"
    echo "  Total chunks: ${NUM_CHUNKS}"
    echo "  Parallel workers: ${PARALLEL_WORKERS}"

    local SLURM_SCRIPT="${RESULTS_DIR}/slurm_${DB_NAME}.sh"

    cat > "${SLURM_SCRIPT}" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=vs_DBNAME
#SBATCH --partition=PARTITION_VAL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=NTASKS_VAL
#SBATCH --mem=MEM_VAL
#SBATCH --time=TIME_VAL
#SBATCH --output=RESULTS_DIR_VAL/logs/DBNAME_%j.out
#SBATCH --error=RESULTS_DIR_VAL/logs/DBNAME_%j.err
ACCOUNT_VAL

module purge
module load cpu
module load gcc
module load openbabel 2>/dev/null

TOTAL_COMPOUNDS=TOTAL_VAL
CHUNK=CHUNK_VAL
WORKERS=WORKERS_VAL
SMI=SMI_VAL
SCREEN_PY=SCREEN_VAL
RECEPTOR_FILE=RECEPTOR_VAL
RES_DIR=RESULTS_DIR_VAL
POSES=POSES_VAL

echo "Starting virtual screening: DBNAME"
echo "Total compounds: ${TOTAL_COMPOUNDS}"
echo "Workers: ${WORKERS}"
echo "Chunks of: ${CHUNK}"
echo ""

# Function to process one chunk
process_chunk() {
    local START=$1
    local CHUNK_ID=$2
    local OUTPUT="${RES_DIR}/DBNAME_chunk_${CHUNK_ID}.tsv"

    if [ -f "$OUTPUT" ]; then
        local LINES=$(wc -l < "$OUTPUT")
        if [ "$LINES" -gt 1 ]; then
            echo "  Chunk ${CHUNK_ID} already done (${LINES} lines), skipping"
            return 0
        fi
    fi

    python3 ${SCREEN_PY} \
        --smi ${SMI} \
        --receptor ${RECEPTOR_FILE} \
        --center CX CY CZ \
        --size SX SY SZ \
        --exhaustiveness EX \
        --start ${START} \
        --count ${CHUNK} \
        --output ${OUTPUT} \
        --save-poses ${POSES} \
        --save-threshold -7.0

    echo "  Chunk ${CHUNK_ID} complete"
}

# Launch chunks in parallel, WORKERS at a time
RUNNING=0
CHUNK_ID=0
for START in $(seq 0 ${CHUNK} $((TOTAL_COMPOUNDS - 1))); do
    process_chunk ${START} ${CHUNK_ID} &
    CHUNK_ID=$((CHUNK_ID + 1))
    RUNNING=$((RUNNING + 1))

    if [ ${RUNNING} -ge ${WORKERS} ]; then
        wait -n 2>/dev/null || wait
        RUNNING=$((RUNNING - 1))
    fi
done

# Wait for all remaining
wait

echo ""
echo "All chunks complete for DBNAME"
echo "Results in: ${RES_DIR}/"
SLURM_EOF

    # Replace placeholders
    sed -i "s|DBNAME|${DB_NAME}|g" "${SLURM_SCRIPT}"
    sed -i "s|PARTITION_VAL|${PARTITION}|g" "${SLURM_SCRIPT}"
    sed -i "s|NTASKS_VAL|${NTASKS}|g" "${SLURM_SCRIPT}"
    sed -i "s|MEM_VAL|${MEM}|g" "${SLURM_SCRIPT}"
    sed -i "s|TIME_VAL|${TIME}|g" "${SLURM_SCRIPT}"
    sed -i "s|ACCOUNT_VAL|#SBATCH ${ACCOUNT}|g" "${SLURM_SCRIPT}"
    sed -i "s|RESULTS_DIR_VAL|${RESULTS_DIR}|g" "${SLURM_SCRIPT}"
    sed -i "s|TOTAL_VAL|${TOTAL}|g" "${SLURM_SCRIPT}"
    sed -i "s|CHUNK_VAL|${CHUNK_SIZE}|g" "${SLURM_SCRIPT}"
    sed -i "s|WORKERS_VAL|${PARALLEL_WORKERS}|g" "${SLURM_SCRIPT}"
    sed -i "s|SMI_VAL|${SMI_FILE}|g" "${SLURM_SCRIPT}"
    sed -i "s|SCREEN_VAL|${SCREEN_SCRIPT}|g" "${SLURM_SCRIPT}"
    sed -i "s|RECEPTOR_VAL|${RECEPTOR}|g" "${SLURM_SCRIPT}"
    sed -i "s|POSES_VAL|${POSES_DIR}|g" "${SLURM_SCRIPT}"
    sed -i "s|CX|${CENTER_X}|g" "${SLURM_SCRIPT}"
    sed -i "s|CY|${CENTER_Y}|g" "${SLURM_SCRIPT}"
    sed -i "s|CZ|${CENTER_Z}|g" "${SLURM_SCRIPT}"
    sed -i "s|SX|${SIZE_X}|g" "${SLURM_SCRIPT}"
    sed -i "s|SY|${SIZE_Y}|g" "${SLURM_SCRIPT}"
    sed -i "s|SZ|${SIZE_Z}|g" "${SLURM_SCRIPT}"
    sed -i "s|EX|${EXHAUSTIVENESS}|g" "${SLURM_SCRIPT}"

    echo "  Submitting: sbatch ${SLURM_SCRIPT}"
    sbatch "${SLURM_SCRIPT}"
    echo ""
}

# ── Main ──

DB_CHOICE="${1:-both}"

if [ ! -f "$RECEPTOR" ]; then
    echo "ERROR: Receptor not found at ${RECEPTOR}"
    echo "Make sure 2ZSH.pdbqt exists in ${PROJECT_DIR}/receptors/"
    exit 1
fi

if [ ! -f "$SCREEN_SCRIPT" ]; then
    echo "ERROR: Screening script not found at ${SCREEN_SCRIPT}"
    exit 1
fi

echo "============================================"
echo "Virtual Screening Submission"
echo "============================================"
echo "Receptor: ${RECEPTOR}"
echo "Box center: ${CENTER_X}, ${CENTER_Y}, ${CENTER_Z}"
echo "Box size: ${SIZE_X} x ${SIZE_Y} x ${SIZE_Z}"
echo "Exhaustiveness: ${EXHAUSTIVENESS}"
echo "Node: ${NTASKS} tasks, ${MEM} RAM, ${TIME}"
echo "Parallel workers: ${PARALLEL_WORKERS}"
echo "============================================"
echo ""

case "$DB_CHOICE" in
    npatlas)
        submit_database "npatlas"
        ;;
    coconut)
        submit_database "coconut"
        ;;
    both)
        submit_database "npatlas"
        submit_database "coconut"
        ;;
    *)
        echo "Usage: bash submit_screening.sh [npatlas|coconut|both]"
        exit 1
        ;;
esac

echo "============================================"
echo "Jobs submitted! Monitor with: squeue -u \$USER"
echo "Check progress: python3 ${PROJECT_DIR}/collect_results.py --progress"
echo "When done: python3 ${PROJECT_DIR}/collect_results.py"
echo "============================================"
