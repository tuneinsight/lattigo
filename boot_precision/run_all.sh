#/bin/sh

go build boot_precision.go

mkdir -p out

./boot_precision -makeplot -paramSet 0 slotcount | tee out/slotcount_p0.tex
./boot_precision -makeplot -paramSet 1 slotcount | tee out/slotcount_p1.tex
./boot_precision -makeplot -paramSet 2 slotcount | tee out/slotcount_p2.tex
./boot_precision -makeplot -paramSet 3 slotcount | tee out/slotcount_p4.tex

./boot_precision -makeplot -logslot 14 -paramSet 0 slotdist | tee out/slotdist_p0_14.tex
./boot_precision -makeplot -logslot 15 -paramSet 0 slotdist | tee out/slotdist_p0_15.tex
./boot_precision -makeplot -logslot 14 -paramSet 1 slotdist | tee out/slotdist_p1_14.tex
./boot_precision -makeplot -logslot 15 -paramSet 1 slotdist | tee out/slotdist_p1_15.tex
./boot_precision -makeplot -logslot 14 -paramSet 2 slotdist | tee out/slotdist_p2_14.tex
./boot_precision -makeplot -logslot 15 -paramSet 2 slotdist | tee out/slotdist_p2_15.tex
./boot_precision -makeplot -logslot 13 -paramSet 3 slotdist | tee out/slotdist_p3_13.tex
./boot_precision -makeplot -logslot 14 -paramSet 3 slotdist | tee out/slotdist_p3_14.tex

./boot_precision -makeplot -logslot 15 -nboot 50 -paramSet 0 successive | tee out/successive_p0.tex
./boot_precision -makeplot -logslot 15 -nboot 50 -paramSet 1 successive | tee out/successive_p1.tex
./boot_precision -makeplot -logslot 15 -nboot 50 -paramSet 2 successive | tee out/successive_p2.tex
./boot_precision -makeplot -logslot 14 -nboot 50 -paramSet 3 successive | tee out/successive_p4.tex