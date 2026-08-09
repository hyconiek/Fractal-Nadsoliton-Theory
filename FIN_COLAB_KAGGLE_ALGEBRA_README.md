# FIN exact-algebra notebook — Colab/Kaggle

Open `FIN_Exact_Algebra_P485_P487_Colab_Kaggle.ipynb` and upload
`FIN_Exact_Algebra_Colab_Kaggle_Bundle.zip` when requested.

The default run is `P486,P485,P487`. The old fourteen-variable P475 route is
disabled by default because it previously ran for more than twelve hours
without returning a basis. Run it only by explicitly changing the program
list to include `P475`.

The bundle contains the repaired coefficient normalization
`sqrt(2)=2-4*alpha**2`, the exact P480 certificate, portable scripts, and a
runner that creates checkpoint ZIP files after every stage. In Colab, enable
Google Drive persistence before extraction if the computation must survive a
browser disconnect. In Kaggle, keep outputs under `/kaggle/working` and
download the final checkpoint before ending the session.

P485 success still requires all five exact ideal remainders to vanish. A
returned P487 relation is not a minimal polynomial until its P473 factor is
isolated, proved irreducible, and independently verified.
