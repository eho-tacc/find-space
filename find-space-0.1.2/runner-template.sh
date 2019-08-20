# Import Agave runtime extensions
. _lib/extend-runtime.sh

# Allow CONTAINER_IMAGE over-ride via local file
if [ -z "${CONTAINER_IMAGE}" ]
then
    if [ -f "./_lib/CONTAINER_IMAGE" ]; then
        CONTAINER_IMAGE=$(cat ./_lib/CONTAINER_IMAGE)
    fi
    if [ -z "${CONTAINER_IMAGE}" ]; then
        echo "CONTAINER_IMAGE was not set via the app or CONTAINER_IMAGE file"
        CONTAINER_IMAGE="sd2e/base:ubuntu17"
    fi
fi

FS_WORK_DIR=${work_dir}
FS_SCAN_DIR=${scan_dir}
FS_OUT_DIR=${out_dir}
FS_PDB_SCREEN_LIST=${pdb_screen_list}
FS_LOG_FP=${log_fp}

# default handlers
if [ -z $FS_WORK_DIR ]; then FS_WORK_DIR="/home/work/"; fi
if [ -z $FS_SCAN_DIR ]; then FS_SCAN_DIR="/home/work/scan/"; fi
if [ -z $FS_OUT_DIR ]; then FS_OUT_DIR="/home/work/out/"; fi
if [ -z $FS_PDB_SCREEN_LIST ]; then FS_PDB_SCREEN_LIST="/home/work/pdb_screen_list.txt"; fi
if [ -z $FS_LOG_FP ]; then FS_LOG_FP="/home/work/log/fs_log.txt"; fi

cmd="python3 /opt/find-space/src/find_space.py --print-params " \
"-w $(FS_WORK_DIR) -s $(FS_SCAN_DIR) -o $(FS_OUT_DIR) " \
"-p $(FS_PDB_SCREEN_LIST) -l $(FS_LOG_FP)"

echo DEBUG=1 container_exec ${CONTAINER_IMAGE} $cmd
DEBUG=1 container_exec ${CONTAINER_IMAGE} $cmd
