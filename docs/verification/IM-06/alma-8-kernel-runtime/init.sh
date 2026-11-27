#!/bin/sh
/bin/busybox mount -t proc proc /proc
/bin/busybox mount -t devtmpfs devtmpfs /dev
/run-probe
printf 'UVLAS_KERNEL_REHEARSAL_COMPLETE\n'
/bin/busybox poweroff -f
