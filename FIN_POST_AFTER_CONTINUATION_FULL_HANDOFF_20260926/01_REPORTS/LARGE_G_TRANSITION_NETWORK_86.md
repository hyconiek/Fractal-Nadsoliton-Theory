# LARGE-G-TRANSITION-NETWORK-86
## Six asymptotic index-one saddle classes organize transitions among twelve one-label minima

Date: 2026-09-26

Status:
- exact consequence of the support classification for the large-g limit;
- branch-to-atlas identifications at g=5 are numerical continuation facts.

The unique size-one D12 support class consists of the twelve label vertices.

For one label i,

    p*=delta_i,

and

    mu_min = A_ii
           ≈ 1.271778833994.

These are the asymptotic index-zero localized minima.

Every two-label support is automatically the equal mixture

    p*=(1/2,1/2).

D12 classifies an unordered pair only by cyclic separation d=1,...,6.

There are therefore exactly six asymptotic index-one saddle classes.

For a distance-d pair,

    mu_d=(A_00+A_0d)/2,

and because both single- and two-label supports have zero 1/g correction,

    V_pair - V_single
      =
      [g/2](mu_min-mu_d)
      -log 2
      + exponentially small terms.

Numerically:

    d=1:
      mu=0.280762177817
      barrier slope=0.495508328088

    d=2:
      mu=0.574156322828
      barrier slope=0.348811255583

    d=3:
      mu=0.721595902727
      barrier slope=0.275091465633

    d=4:
      mu=0.709506943817
      barrier slope=0.281135945088

    d=5:
      mu=0.612537415303
      barrier slope=0.329620709345

    d=6:
      mu=0.561776644984
      barrier slope=0.355001094505.

Thus the cheapest large-g index-one saddle is the distance-3 pair.

Numerical continuation identifies the g=5 atlas index-one branch i=1 with
this d=3 asymptotic support.

The next cheapest class is d=4; the g=5 secondary index-one branch i=2
continues to a d=4 pair.

The reflection index-one branch i=4 continues to a d=5 pair.

The d=1,d=2,d=6 index-one branches are asymptotically guaranteed but are not
present among the known g=5 index-one atlas orbits, so they must be created or
connected through events above g=5.

Graphically, d=3 edges alone preserve label residue modulo 3 and split the
twelve minima into three four-cycles.

Adding the next-cheapest d=4 edges makes the Cayley graph on Z12 connected
because gcd(3,4,12)=1.

Therefore the asymptotic communication hierarchy is:

    first:  low barriers organize four-state cycles through d=3 moves;

    second: d=4 saddles connect those sectors into one global twelve-minimum
            network.

This is a concrete metastable graph prediction of the large-g FIN landscape.
