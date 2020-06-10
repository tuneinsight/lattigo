#/bin/sh

go build boot_precision.go

mkdir -p out

./boot_precision -makeplot -logslot 4 -nboot 3 -paramSet 1 slotcount | tee out/slotcount_p1.tex

./boot_precision -makeplot -logslot 4 -nboot 3 -paramSet 1 successive | tee out/successive_p1.tex

./boot_precision -makeplot -logslot 4 -nboot 3 -paramSet 1 slotdist | tee out/slotdist_p1.tex