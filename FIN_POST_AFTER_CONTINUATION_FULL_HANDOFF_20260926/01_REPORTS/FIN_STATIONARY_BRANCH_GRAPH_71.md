# FIN-STATIONARY-BRANCH-GRAPH-71
## A first connected branch graph joining fold localization, D3 splitting, and the g=5 atlas

Date: 2026-09-26

Status:
- synthesis of certified local events and numerical continuation;
- graph edges marked by proof level;
- no claim of exhaustive full-X7 stationary topology.

## 1. Certified nodes

### F_low
Simple fold:

    g ≈ 3.51564471684

Proof level:
    interval-assisted controlled local theorem (R7P-031 / MP7-026).

It connects:
    index-0 localized minimum
    and
    index-1 saddle.

### C_D3
D3 transverse crossing:

    g ≈ 5.17184183194

Proof level:
    interval-assisted linear crossing (MP7-039)
    plus interval-separated nonzero cubic coefficient (report 66).

Base index changes locally:

    index 3 -> index 1.

Two D3 daughter triplets have index 2 at birth.

## 2. New numerical nodes

### F_upper

    g ≈ 5.172231474684

connects:
    index-2 D3 daughter sheet
    and
    index-1 main saddle sheet.

### F_second

    g ≈ 4.395526393543

connects:
    index-2 second-daughter sheet
    and
    index-1 secondary saddle sheet.

These folds are currently numerical.

## 3. g=5 orbit identifications

The connected branches account for at least the following R7P-037 atlas
orbits:

    orbit i=0:
      index 0,
      stable localized branch;

    orbit i=7:
      index 1,
      main saddle branch;

    orbit i=9:
      index 2,
      second D3 daughter pre-fold;

    orbit i=2:
      index 1,
      second D3 daughter post-fold.

All four matches are exact up to D12 action at numerical precision
(<2e-13 in feature coordinates).

## 4. Branch graph

A compact schematic is:

                               C_D3
                         g≈5.171841832
                         /             \
                        /               \
              triplet +                 triplet -
               idx 2                     idx 2
                 |                         |
          F_upper≈5.17223147        F_second≈4.39552639
                 |                         |
               idx 1                     idx 1
                 |                         |
            atlas i=7                 atlas i=2
                 |
          F_low≈3.51564472
              /     \
           idx 1    idx 0
                    |
                 atlas i=0

The second triplet's pre-fold sheet is atlas i=9 at g=5.

## 5. Conceptual result

The FIN stationary landscape is beginning to look less like a collection of
unrelated numerical roots and more like a symmetry-controlled branch network.

In this network:

- simple folds create/annihilate pairs of different Morse index;
- D3 symmetry crossings redistribute index through a two-dimensional critical
  representation;
- D12 orbit/stabilizer changes control multiplicity;
- several apparently unrelated g=5 roots are connected pieces of this same
  network.

## 6. Remaining global obligation

This graph is not exhaustive.

The g=5 atlas contains additional index-0,1,2,3 orbits not yet assigned to
validated connected components.

The next useful campaign is therefore not another random root census but:

    BRANCH-GRAPH-COMPLETION-72

Continue the remaining g=5 orbit representatives in g, detect their folds and
symmetry crossings, and build a quotient graph of stationary D12 orbits.

For every edge record:
- orbit stabilizer;
- Morse index;
- energy;
- branch endpoints/crossings;
- proof level.

The target is a finite stationary-orbit graph, not merely a point census.
