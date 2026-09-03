#!/usr/bin/env python3
"""Combinatorial product-refinement helpers for FIN cell-complex studies."""
from __future__ import annotations

import itertools


def gf2_rank(columns: list[int]) -> int:
    pivots={}
    for value in columns:
        while value:
            pivot=value.bit_length()-1
            if pivot in pivots:value^=pivots[pivot]
            else:pivots[pivot]=value;break
    return len(pivots)


def base_complex():
    vertices=list(range(12))
    edges=list(itertools.combinations(vertices,2))
    faces=[tuple(t) for t in itertools.combinations(vertices,3)]
    return vertices,edges,faces


def refine(vertices,edges,faces):
    index={(v,b):2*i+b for i,v in enumerate(vertices) for b in (0,1)}
    new_vertices=list(range(2*len(vertices)))
    horizontal=[]
    for a,b in edges:
        for fiber in (0,1):horizontal.append(tuple(sorted((index[(a,fiber)],index[(b,fiber)]))))
    vertical=[tuple(sorted((index[(v,0)],index[(v,1)]))) for v in vertices]
    new_edges=sorted(horizontal+vertical)
    lifted=[]
    for face in faces:
        for fiber in (0,1):lifted.append(tuple(index[(v,fiber)] for v in face))
    squares=[]
    for a,b in edges:squares.append((index[(a,0)],index[(b,0)],index[(b,1)],index[(a,1)]))
    return new_vertices,new_edges,lifted+squares,{"lifted_faces":len(lifted),"new_squares":len(squares)}


def boundary_bits(edges,faces):
    edge_index={tuple(sorted(e)):i for i,e in enumerate(edges)}
    output=[]
    for face in faces:
        value=0
        for a,b in zip(face,face[1:]+face[:1]):value^=1<<edge_index[tuple(sorted((a,b)))]
        output.append(value)
    return output

