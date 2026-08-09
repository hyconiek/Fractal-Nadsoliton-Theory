# FIN exact algebra with Singular

Use `FIN_Exact_Algebra_Singular_Colab_Kaggle.ipynb` with
`FIN_Exact_Algebra_Singular_Bundle.zip`.

The notebook keeps a normal Python kernel. SymPy constructs and audits the
exact rational input; **Singular** performs P485 ideal reduction and P487
elimination. The recommended campaign is:

`P486 -> PREPARE -> P485(Singular) -> P487(Singular)`

P475 is available only as a separate explicit action. It is not part of the
recommended run.

If `Singular` is absent, the notebook contains an explicit installation cell
for Colab/Kaggle. A SageMath environment can also run the same notebook
because Sage normally includes Singular.

Every stage writes status JSON, logs, exact Singular input, and checkpoint
ZIP files. P485 passes only if all five printed normal forms are exactly zero
and P486 premises pass. P487 never labels a returned relation a minimal
polynomial without factor isolation, irreducibility, and exact verification.
