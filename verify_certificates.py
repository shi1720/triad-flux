"""Dependency-free checker for saved rational phase-compatibility certificates.

This validates inequalities for the saved cosine polynomial. Independent FFT
and convolution tests validate the fluid-to-polynomial conversion separately.
"""
import json
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path

def verify(item):
    rows=item['rows'];coeff=list(map(Q,item['coefficients']));used=[Q(0)]*len(rows)
    assert len(rows)==len(coeff) and all(coeff)
    tri=sum(map(abs,coeff),Q(0));assert tri==Q(item['triangle'])
    packed=Q(0);general=Q(0)
    for block in item['cycles']+item.get('general_relations',[]):
        rel=block['relation'];m=Q(block['mass']);assert m>0
        assert len(set(i for i,_ in rel))==len(rel)
        assert all(0<=i<len(rows) and isinstance(z,int) and z for i,z in rel)
        assert all(sum(z*rows[i][j] for i,z in rel)==0 for j in range(len(rows[0])))
        assert sum(z for i,z in rel if coeff[i]<0)%2==1
        for i,_ in rel:used[i]+=m
        if 'penalty' in block:
            p=Q(2,sum(z*z for _,z in rel));assert p==Q(block['penalty']);general+=m*p
        else:
            assert len(rel)==4 and all(abs(z)==1 for _,z in rel);packed+=m
    assert all(u<=abs(c) for u,c in zip(used,coeff))
    assert packed==Q(item['packed_mass']) and general==Q(item['general_deficit'])
    # sqrt(2) is bounded without floating-point arithmetic.
    scale=10**30;root=isqrt(2*scale*scale);lo=Q(root,scale);hi=Q(root+1,scale)
    assert lo*lo<=2<=hi*hi
    upper=tri-4*packed+2*hi*packed-general
    return upper

if __name__=='__main__':
    items=json.loads((Path(__file__).parent/'results/packing_certificates.json').read_text())
    for item in items:
        u=verify(item);cond=item['condition']
        print(f"L={cond['L']}, p={cond['decay']}, {cond['observable']}: VERIFIED upper <= {float(u):.14g}")
    print(f'{len(items)} exact rational certificates verified.')
