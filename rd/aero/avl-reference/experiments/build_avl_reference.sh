#!/bin/sh
set -eu
uvlas_root="$PWD"
uvlas_avl="$uvlas_root/build/reference-avl"
uvlas_prefix="$uvlas_avl/toolchain"
export PATH="$uvlas_prefix/bin:$PATH"
export LD_LIBRARY_PATH="$uvlas_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
uvlas_fc="$uvlas_prefix/bin/x86_64-conda-linux-gnu-gfortran"
uvlas_cc="$uvlas_prefix/bin/x86_64-conda-linux-gnu-gcc"
"$uvlas_fc" --version > "$uvlas_avl/compiler-version.txt"
"$uvlas_cc" --version >> "$uvlas_avl/compiler-version.txt"
"$uvlas_root/build/micromamba-bootstrap/bin/micromamba" --no-rc -r "$uvlas_avl/mamba-root" list -p "$uvlas_prefix" --explicit > "$uvlas_avl/toolchain-explicit.txt"
uvlas_original="$uvlas_avl/source-baseline"
make -C "$uvlas_original/plotlib" -f Makefile.all lib -j2 \
 FC="$uvlas_fc" CC="$uvlas_cc" PLTLIB=libPlt_gDP.a \
 FFLAGS='-O2 -fdefault-real-8 -fallow-argument-mismatch' \
 CFLAGS="-O2 -DUNDERSCORE -DDBL_ARGS -I$uvlas_prefix/include" \
 CFLAGS0="-O0 -DUNDERSCORE -DDBL_ARGS -I$uvlas_prefix/include" \
 INCDIR="-I$uvlas_prefix/include" > "$uvlas_avl/build-plotlib.log" 2>&1
make -C "$uvlas_original/eispack" -f Makefile.gfortran -j2 \
 FC="$uvlas_fc" FLG='-O2 -fallow-argument-mismatch' DP=-fdefault-real-8 \
 > "$uvlas_avl/build-eispack.log" 2>&1
for uvlas_variant in baseline measured core0; do
 uvlas_dir="$uvlas_avl/source-$uvlas_variant"
 uvlas_extra=""
 if [ "$uvlas_variant" != baseline ]; then
  "$uvlas_fc" -c -O -fdefault-real-8 -fallow-argument-mismatch \
   -I"$uvlas_dir/src" -I"$uvlas_dir/bin" "$uvlas_dir/src/avl_heap_inc.f90" \
   -J"$uvlas_dir/bin" -o "$uvlas_dir/bin/avl_heap_inc.o"
  "$uvlas_fc" -c -O -fdefault-real-8 -fallow-argument-mismatch \
   -I"$uvlas_dir/src" -I"$uvlas_dir/bin" "$uvlas_dir/src/uvl_measurement.f" \
   -o "$uvlas_dir/bin/uvl_measurement.o"
  uvlas_extra="uvl_measurement.o"
 fi
 make -C "$uvlas_dir/bin" -f Makefile.gfortranDP avl \
  FC="$uvlas_fc" PLTOBJ="$uvlas_original/plotlib/libPlt_gDP.a" \
  EIGOBJ="$uvlas_original/eispack/libeispack.a" ALIBS="$uvlas_extra" \
  PLTLIB="-L$uvlas_prefix/lib -lX11 -Wl,-rpath,$uvlas_prefix/lib" \
  > "$uvlas_avl/build-$uvlas_variant.log" 2>&1
done
