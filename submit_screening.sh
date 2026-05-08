#!/bin/bash
#=====================================================================
# Virtual Screening Submission Script for SDSC Expanse
# Screens NP Atlas and/or COCONUT databases using AutoDock Vina
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

# Chunk size: compounds per Slurm task
CHUNK_SIZE=200

# Slurm settings
PARTITION="shared"
TIME="04:00:00"
CPUS=1
MEM="4G"
ACCOUNT=""  # Set your allocation account if required, e.g., ACCOUNT="--account=TG-CHE123456"

mkdir -p "${RESULTS_DIR}" "${POSES_DIR}"

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
    local NUM_TASKS=$(( (TOTAL + CHUNK_SIZE - 1) / CHUNK_SIZE ))

    echo "Database: ${DB_NAME}"
    echo "  Compounds: ${TOTAL}"
    echo "  Chunk size: ${CHUNK_SIZE}"
    echo "  Array tasks: ${NUM_TASKS}"

    local SLURM_SCRIPT="${RESULTS_DIR}/slurm_${DB_NAME}.sh"

    cat > "${SLURM_SCRIPT}" << SLURM_EOF
#!/bin/bash
#SBATCH --job-name=vs_${DB_NAME}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --array=0-$((NUM_TASKS - 1))%100
#SBATCH --output=${RESULTS_DIR}/logs/${DB_NAME}_%a.out
#SBATCH --error=${RESULTS_DIR}/logs/${DB_NAME}_%a.err
${ACCOUNT}

module purge
module load cpu
module load gcc
module load openbabel 2>/dev/null

START=\$(( SLURM_ARRAY_TASK_ID * ${CHUNK_SIZE} ))
OUTPUT="${RESULTS_DIR}/${DB_NAME}_chunk_\${SLURM_ARRAY_TASK_ID}.tsv"

echo "Task \${SLURM_ARRAY_TASK_ID}: compounds \${START} to \$(( START + ${CHUNK_SIZE} - 1 ))"

python3 ${SCREEN_SCRIPT} \\
    --smi ${SMI_FILE} \\
    --receptor ${RECEPTOR} \\
    --center ${CENTER_X} ${CENTER_Y} ${CENTER_Z} \\
    --size ${SIZE_X} ${SIZE_Y} ${SIZE_Z} \\
    --exhaustiveness ${EXHAUSTIVENESS} \\
    --start \${START} \\
    --count ${CHUNK_SIZE} \\
    --output \${OUTPUT} \\
    --save-poses ${POSES_DIR} \\
    --save-threshold -7.0

echo "Task \${SLURM_ARRAY_TASK_ID} complete"
SLURM_EOF

    mkdir -p "${RESULTS_DIR}/logs"
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
    echo "Download it: curl -sL https://raw.githubusercontent.com/anirudhatalmale6-alt/autogrow-filter-fix/main/screen_database.py -o ${SCREEN_SCRIPT}"
    exit 1
fi

echo "============================================"
echo "Virtual Screening Submission"
echo "============================================"
echo "Receptor: ${RECEPTOR}"
echo "Box center: ${CENTER_X}, ${CENTER_Y}, ${CENTER_Z}"
echo "Box size: ${SIZE_X} x ${SIZE_Y} x ${SIZE_Z}"
echo "Exhaustiveness: ${EXHAUSTIVENESS}"
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
echo "When done, run: python3 ${PROJECT_DIR}/collect_results.py"
echo "============================================"
