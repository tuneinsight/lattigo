#/bin/sh

go build boot_precision.go

mkdir -p out

#./boot_precision -makeplot -logslot 15 -nboot 1 -paramSet 0 slotcount | tee out/slotcount_p0.tex
#./boot_precision -makeplot -logslot 15 -nboot 1 -paramSet 1 slotcount | tee out/slotcount_p1.tex
#./boot_precision -makeplot -logslot 15 -nboot 1 -paramSet 2 slotcount | tee out/slotcount_p2.tex
#./boot_precision -makeplot -logslot 15 -nboot 1 -paramSet 3 slotcount | tee out/slotcount_p3.tex
#./boot_precision -makeplot -logslot 15 -nboot 1 -paramSet 4 slotcount | tee out/slotcount_p4.tex
./boot_precision -makeplot -logslot 15 -nboot 1 -paramSet 5 slotcount | tee out/slotcount_p5.tex
./boot_precision -makeplot -logslot 15 -nboot 1 -paramSet 6 slotcount | tee out/slotcount_p6.tex