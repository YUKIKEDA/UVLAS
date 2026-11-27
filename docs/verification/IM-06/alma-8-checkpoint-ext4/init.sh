#!/bin/sh
set -eu
export PATH=/bin LC_ALL=C HOME=/tmp
b=/bin/busybox

halt() {
  result=$?
  if [ "$result" -ne 0 ]; then
    printf 'UVLAS_CHECKPOINT_GUEST_FAILED=%s\n' "$result"
  fi
  "$b" poweroff -f
}
trap halt EXIT
"$b" mount -t proc proc /proc
"$b" mount -t devtmpfs devtmpfs /dev
"$b" uname -a
for module in virtio_blk mbcache jbd2 libcrc32c ext4; do
  "$b" insmod "/modules/$module.ko"
done
"$b" mount -t ext4 /dev/vda /save
"$b" grep -F '/dev/vda /save ext4 rw,' /proc/mounts
"$b" sha256sum /tests/* /lib64/libc.so.6 /lib64/libstdc++.so.6
cd /save
"$b" mkdir -p /save/tmp
export TMPDIR=/save/tmp

printf '=== existing directory tests at %s ===\n' "$PWD"
/tests/checkpoint_directory
printf '=== existing publication tests at TMPDIR=%s ===\n' "$TMPDIR"
/tests/checkpoint_publication
printf '=== interrupted publication and recovery under the same parent ===\n'
/tests/recovery_probe /save/recovery-proof
"$b" find /save/recovery-proof -type f | "$b" sort
cd /
"$b" sync
"$b" umount /save
printf 'UVLAS_CHECKPOINT_GUEST_COMPLETE\n'
