#/bin/sh

go build boot_precision.go

mkdir -p out

./boot_precision -makeplot -logslot 14 -paramSet 0 slotdist | tee out/slotdist_p0_14.tex
./boot_precision -makeplot -logslot 15 -paramSet 0 slotdist | tee out/slotdist_p0_15.tex

./boot_precision -makeplot -logslot 14 -paramSet 1 slotdist | tee out/slotdist_p1_14.tex
./boot_precision -makeplot -logslot 15 -paramSet 1 slotdist | tee out/slotdist_p1_15.tex

./boot_precision -makeplot -logslot 14 -paramSet 2 slotdist | tee out/slotdist_p2_14.tex
./boot_precision -makeplot -logslot 15 -paramSet 2 slotdist | tee out/slotdist_p2_15.tex

./boot_precision -makeplot -logslot 14 -paramSet 3 slotdist | tee out/slotdist_p3_14.tex
./boot_precision -makeplot -logslot 15 -paramSet 3 slotdist | tee out/slotdist_p3_15.tex

./boot_precision -makeplot -logslot 13 -paramSet 4 slotdist | tee out/slotdist_p4_13.tex
./boot_precision -makeplot -logslot 14 -paramSet 4 slotdist | tee out/slotdist_p4_14.tex