#!/usr/bin/env python3
"""FIN ST2352--ST2366: exact regular triangular-prism audit."""
import itertools
import json
import math

import sympy as sp

from fin_st2172_st2261_common import d12_permutations, strict_A_W, write_packet, write_round


LO, HI = 2352, 2366
NAMES = [
    "PrismDefinition", "ExactCoordinates", "ProductCellDecomposition",
    "OrientedIncidence", "ChainIdentity", "BoundaryHomology", "SolidHomology",
    "GraphSpectrum", "UnweightedHodgeSpectrum", "WeightedHodgeSpectrum",
    "SymmetryAndSelector", "FINRefinementRole", "TernaryTransportRole",
    "PhysicalScaleBoundary", "FinalPrismVerdict",
]
EDGES = [(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3),(1,4),(2,5)]
FACES = [(0,2,1),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5)]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k-LO], status, boundary, payload)


def matrices():
    ei={e:i for i,e in enumerate(EDGES)}
    d0=sp.zeros(9,6); d1=sp.zeros(5,9)
    for r,(i,j) in enumerate(EDGES): d0[r,i],d0[r,j]=-1,1
    for r,face in enumerate(FACES):
        for a,b in zip(face,face[1:]+face[:1]):
            d1[r,ei[(min(a,b),max(a,b))]] += 1 if a<b else -1
    d2=sp.ones(1,5)
    return d0,d1,d2


def eigdict(M):
    return {str(k):int(v) for k,v in M.eigenvals().items()}


def int_matrix(M):
    return [[int(v) for v in row] for row in M.tolist()]


def embedding_scan():
    _,W=strict_A_W(); seen=set(); rows=[]
    for six_tuple in itertools.combinations(range(12),6):
        six=set(six_tuple); root=min(six)
        for pair in itertools.combinations(sorted(six-{root}),2):
            a=tuple(sorted((root,)+pair)); b=tuple(sorted(six-set(a)))
            base=set(itertools.combinations(a,2))|set(itertools.combinations(b,2))
            for matching in itertools.permutations(b):
                edges=base|{tuple(sorted((a[i],matching[i]))) for i in range(3)}
                key=tuple(sorted(edges))
                if key in seen: continue
                seen.add(key); vals=[W[i,j] for i,j in key]
                rows.append((sum(vals),sum(math.log(v) for v in vals),key))
    output={"embeddings":len(rows)}; best_sets=[]
    for index,name in ((0,"sum"),(1,"log_product")):
        best=max(row[index] for row in rows)
        tops=[row for row in rows if abs(row[index]-best)<1e-13]
        top_keys={row[2] for row in tops}; best_sets.append(top_keys)
        orbit=set()
        for _,_,p in d12_permutations():
            orbit.add(tuple(sorted(tuple(sorted((p[i],p[j]))) for i,j in tops[0][2])))
        levels=sorted(set(round(row[index],14) for row in rows),reverse=True)
        output[name]={"best":best,"degeneracy":len(tops),"is_one_D12_orbit":top_keys==orbit,
                      "gap_to_next_level":levels[0]-levels[1],"representative_edges":[list(e) for e in tops[0][2]]}
    output["same_optimum_for_sum_and_product"]=best_sets[0]==best_sets[1]
    return output


def main():
    d0,d1,d2=matrices()
    scan=embedding_scan()
    adjacency=sp.zeros(6,6)
    for i,j in EDGES: adjacency[i,j]=adjacency[j,i]=1
    graph_L=3*sp.eye(6)-adjacency
    L1=d0*d0.T+d1.T*d1
    x={}
    x["ST2352"]=packet(2352,"Defined","A right regular triangular prism has equilateral bases; height h remains free. The uniform prism additionally sets h=a.",{
        "vertices":6,"edges":9,"faces":5,"triangular_faces":2,"rectangular_faces":3,"volume_cells":1})
    x["ST2353"]=packet(2353,"Proven geometry","Coordinates use side a and height h; h=a makes all nine edges equal and side faces square.",{
        "coordinates":["(0,0,-h/2)","(a,0,-h/2)","(a/2,sqrt(3)a/2,-h/2)","(0,0,h/2)","(a,0,h/2)","(a/2,sqrt(3)a/2,h/2)"],
        "volume":"sqrt(3)*a^2*h/4","surface_area":"sqrt(3)*a^2/2+3*a*h"})
    x["ST2354"]=packet(2354,"Proven product structure","Exact cell counts of filled triangle times interval.",{
        "identity":"Delta2 x I","C0":6,"C1":"3*2+3*1=9","C2":"1*2+3*1=5","C3":1,
        "faces":["2 horizontal triangles","3 vertical squares"]})
    x["ST2355"]=packet(2355,"Constructed exact", "Integral primitive incidence with outward face orientations.",{
        "d0":int_matrix(d0),"d1":int_matrix(d1),"d2":int_matrix(d2),
        "rank_d0":int(d0.rank()),"rank_d1":int(d1.rank()),"rank_d2":int(d2.rank())})
    x["ST2356"]=packet(2356,"Proven exact","Boundary of a boundary vanishes.",{
        "d1_d0_zero":d1*d0==sp.zeros(5,6),"d2_d1_zero":d2*d1==sp.zeros(1,9)})
    x["ST2357"]=packet(2357,"Proven exact","The five-face boundary complex is a cellulation of S2.",{
        "boundary_betti":[1,0,1],"Euler_characteristic":6-9+5})
    x["ST2358"]=packet(2358,"Proven exact","Adding the one prism volume fills S2 and gives a 3-ball.",{
        "solid_betti":[1,0,0,0],"Euler_characteristic":6-9+5-1})
    x["ST2359"]=packet(2359,"Proven exact","Prism graph is C3 Cartesian-product K2; eigenvalues are additive under the product.",{
        "adjacency_spectrum":eigdict(adjacency),"graph_Laplacian_spectrum":eigdict(graph_L),
        "graph_Laplacian_factor":"lambda*(lambda-2)*(lambda-3)^2*(lambda-5)^2"})
    x["ST2360"]=packet(2360,"Proven exact","Unweighted cellular one-form Hodge Laplacian.",{
        "L1_spectrum":eigdict(L1),"characteristic_polynomial":"(lambda-2)*(lambda-3)^3*(lambda-5)^5","zero_modes":0})
    x["ST2361"]=packet(2361,"Proven exact","D3h-natural face weights leave two positive moduli: t on triangles and q on squares.",{
        "gradient_sector":[2,3,3,5,5],"coexact_sector":["3*t","5*q (multiplicity 2)","3*t+2*q"],
        "unweighted_t_q":[1,1],"unique_ratio_selected":False})
    x["ST2362"]=packet(2362,"Proven no-go","D3h has order 12 and includes base exchange/reflection. It preserves a two-orbit face split and cannot choose top over bottom.",{
        "symmetry_group":"D3h ~= D3 x C2","order":12,"face_orbits":2,"canonical_layer_arrow":False,"QW2191_discharged":False})
    x["ST2363"]=packet(2363,"Conditional positive / exhaustive no-selector","For every selected base triangle, product incidence supplies lifted cells; exhaustive strict-W scoring does not choose a unique global prism.",{
        "fixes_vertical_square_incidence":True,"selects_base_triangle":False,"strict_global_cell_complex_selected":False,
        "embedding_scan":scan})
    x["ST2364"]=packet(2364,"Conditional positive","Triangle data may be transported between layers; square faces carry edge-pair transport. No new connected third cumulant is created.",{
        "transports_existing_tau":True,"generates_irreducible_ternary_source":False,"top_minus_bottom_requires_orientation":True})
    x["ST2365"]=packet(2365,"Boundary","Geometry retains aspect ratio h/a and overall scale a; the cellular Hodge model retains t/q.",{
        "dimensionless_shape_modulus":"h/a unless uniform h=a","absolute_scale":"a","Hodge_ratio":"t/q","physical_values_derived":False})
    x["ST2366"]=packet(2366,"Final verdict","The triangular prism is the correct local refinement cell and repairs conditional product-cell bookkeeping, but it does not source the base face, layer arrow, ternary law, Hodge ratio, scale or physics.",{
        "recommended_role":"canonical conditional local cell Delta2 x I","strict_missing_objects_closed":[],
        "next_test":"audit whether the 12-state carrier admits a strict prism decomposition or only many symmetry-related decompositions","strict_ToE_closure":False})
    write_round(LO,HI,x)


if __name__=="__main__":main()
