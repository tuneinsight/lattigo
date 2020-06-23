#/bin/sh

go build boot_precision.go

mkdir -p out

#./boot_precision -makeplot -logslot 10 -nboot 1 -paramSet 0 slotdist | tee out/slotdist_p0.tex
#./boot_precision -makeplot -logslot 10 -nboot 1 -paramSet 1 slotdist | tee out/slotdist_p1.tex
#./boot_precision -makeplot -logslot 10 -nboot 1 -paramSet 2 slotdist | tee out/slotdist_p2.tex
#./boot_precision -makeplot -logslot 10 -nboot 1 -paramSet 3 slotdist | tee out/slotdist_p3.tex
./boot_precision -makeplot -logslot 10 -nboot 1 -paramSet 4 slotdist | tee out/slotdist_p4.tex