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

# mkdirs if they do not exist
mkdir -p ${work_dir} ${scan_dir} ${out_dir} ${pdb_screen_list} "$(dirname ${log_fp})"

cmd="python3 /opt/find-space/src/find_space.py --print-params -w ${work_dir} -s ${scan_dir} -o ${out_dir} -p ${pdb_screen_list} -l ${log_fp}"
# cmd="ls $FS_WORK_DIR"

echo DEBUG=1 container_exec ${CONTAINER_IMAGE} $cmd
DEBUG=1 container_exec ${CONTAINER_IMAGE} $cmd
