# Z3-VS-Z4-QUOTIENT-COMPETITION-220
## Both CRT factors define exact effective quotients; the Z4 relaxation gap is slightly smaller for N=3..8

Date: 2026-09-26

Status:
exact quotient algebra for the 12-state circulant effective chain;
finite-N numerical comparison for N=3..8.

Because

    Z12 ≅ Z3×Z4,

the circulant localized-state generator is strongly lumpable onto BOTH:
- j mod 3;
- j mod 4.

The Z3 quotient has one nonzero eigenvalue pair k=4,8.

The Z4 quotient has modes k=3,6,9.

Measured gaps:


    N=3:
      gap_Z3=0.132005956666
      gap_Z4=0.127456889351
      gap_Z4/gap_Z3=0.965539

    N=4:
      gap_Z3=0.0718543830295
      gap_Z4=0.068772570252
      gap_Z4/gap_Z3=0.957110

    N=5:
      gap_Z3=0.0402247566571
      gap_Z4=0.0382058190385
      gap_Z4/gap_Z3=0.949809

    N=6:
      gap_Z3=0.0226189121545
      gap_Z4=0.0213311114407
      gap_Z4/gap_Z3=0.943065

    N=7:
      gap_Z3=0.0126419110629
      gap_Z4=0.0118419419585
      gap_Z4/gap_Z3=0.936721

    N=8:
      gap_Z3=0.00698992202668
      gap_Z4=0.00650614914865
      gap_Z4/gap_Z3=0.930790

Over the entire tested range the Z4 quotient has the slightly smaller global
spectral gap.

Therefore the earlier phrase

    "fast Z4 fiber -> slow Z3 base"

is NOT established at N<=8.

The barrier ordering d3<d4 alone is insufficient to rank quotient relaxation
times because:
- several shell channels contribute;
- C3 and C4 have different spectral geometry;
- memory renormalizes each shell differently.

This is an important correction.
