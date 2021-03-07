#/bin/sh

go build boot_precision.go

mkdir -p out

./boot_precision -makeplot -logslot 15 -nboot 50 -paramSet 0 successive | tee out/successive_p0.tex
./boot_precision -makeplot -logslot 15 -nboot 50 -paramSet 1 successive | tee out/successive_p1.tex
./boot_precision -makeplot -logslot 15 -nboot 50 -paramSet 2 successive | tee out/successive_p2.tex
./boot_precision -makeplot -logslot 15 -nboot 50 -paramSet 3 successive | tee out/successive_p3.tex
./boot_precision -makeplot -logslot 14 -nboot 50 -paramSet 4 successive | tee out/successive_p4.tex