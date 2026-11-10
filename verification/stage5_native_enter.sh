#!/bin/bash
# Enter the retained CI root filesystem with the D-drive study bound at /stage5.
set -eu
old=/mnt/d/DivertorOptimization/studies/firefly-baseline-20261105
study=$1
code=$2
shift 2
marker="$study/native-physics/active-library.txt"
if [ -f "$marker" ] && [ "${STAGE5_PRIVATE_NS:-0}" != 1 ]; then
    exec unshare --mount --propagation private env STAGE5_PRIVATE_NS=1 bash "$0" "$study" "$code" "$@"
fi
mnt=/mnt/firefly20261105
r="$mnt/bundle/rootfs"
exec 9>/run/lock/nf-stage5-runtime.lock
flock 9
mkdir -p "$mnt"
mountpoint -q "$mnt" || mount -o loop "$old/runtime.ext4" "$mnt"
mkdir -p "$r/work" "$r/stage5" "$r/stage5code"
mountpoint -q "$r/work" || mount --bind "$old" "$r/work"
libdir=/work/install/lib:/usr/local/lib
thermal=/work/install/lib/libFIREFLY.so
if [ -f "$marker" ]; then
    test "${STAGE5_PRIVATE_NS:-0}" = 1
    mount --bind "$study" "$r/stage5"
    mount --bind "$code" "$r/stage5code"
    relative=$(cat "$marker")
    [[ "$relative" =~ ^native-physics/[A-Za-z0-9-]+/libFIREFLY\.so$ ]]
    thermal="/stage5/$relative"
    test -f "$r$thermal"
    libdir="${thermal%/*}:$libdir"
else
    mountpoint -q "$r/stage5" || mount --bind "$study" "$r/stage5"
    mountpoint -q "$r/stage5code" || mount --bind "$code" "$r/stage5code"
fi
test "$study" -ef "$r/stage5"
test "$code" -ef "$r/stage5code"
for fs in dev proc sys; do
    mountpoint -q "$r/$fs" || mount --bind "/$fs" "$r/$fs"
done
flock -u 9
exec 9>&-
exec chroot "$r" /usr/bin/env PATH=/usr/local/bin:/usr/bin:/bin LD_LIBRARY_PATH="$libdir" STAGE5_THERMAL_LIBRARY="$thermal" PYTHONPATH=/usr/local/lib/python:/work/install/lib/python PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 "$@"
