# Replay

From the repository root:

    OPENBLAS_NUM_THREADS=1 python fin_pde_parameter_campaign_20261108/campaign.py

Requires NumPy and SymPy. The frozen canonical reconstruction is imported
from the previous pilot. Thirty complete nonlinear runs are computed.
results.json retains every twentieth recorded diagnostic and each final
point, plus full-history summary statistics. Runtime is several minutes.
The time step is reduced to reach T=epsilon^-2 exactly. Parameters and the
last-three-quarters inverse-action slope estimator are unchanged from
the pilot. This is the same author's independent numerical route relative
to the reduced action equation, not independent-team replication.
