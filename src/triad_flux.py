"""Exact Fourier flux coefficients and constructive four-triad certificates.

Domain (R / 2pi Z)^2. Fourier coefficients of vorticity. Mean spatial
integral normalized to one. Positive flux is loss from |k|^2 <= cutoff2.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd, isqrt
import numpy as np

def neg(p): return tuple(-x for x in p)
def add(p,q): return tuple(x+y for x,y in zip(p,q))
def norm2(p): return sum(x*x for x in p)
def cross(p,q): return p[0]*q[1]-p[1]*q[0]
def canonical(p): return p if p[0]>0 or (p[0]==0 and p[1]>0) else neg(p)
def phase_row(modes, terms):
    row=[0]*len(modes)
    for m,s in terms:
        cm=canonical(m);row[modes.index(cm)]+=s*(1 if m==cm else -1)
    if not any(row): return None
    if next(x for x in row if x)<0: row=[-x for x in row]
    return tuple(row)

def coefficients(modes, cutoff2, observable='energy', amplitudes=None):
    """Rational c_e in Pi=sum c_e cos(M_e phi), including conjugate modes.

    amplitudes may be exact rational values; omitted means unit vorticity.
    All interactions among occupied modes are enumerated, not a handpicked ODE.
    """
    modes=list(map(tuple,modes))
    if any(p!=(canonical(p)) or p==(0,0) for p in modes) or len(set(modes))!=len(modes):
        raise ValueError('distinct nonzero canonical half-plane modes required')
    if observable not in ('energy','enstrophy'):raise ValueError('unknown observable')
    if amplitudes is not None and len(amplitudes)!=len(modes): raise ValueError('one amplitude per mode required')
    amps={p:Q(a) for p,a in zip(modes,amplitudes if amplitudes is not None else [1]*len(modes))}
    if any(a<=0 for a in amps.values()): raise ValueError('positive amplitudes required')
    signed=set(modes)|{neg(p) for p in modes};out={}
    for k in sorted(signed):
        if norm2(k)>cutoff2: continue
        wk=Q(1,norm2(k)) if observable=='energy' else Q(1)
        for p in sorted(signed):
            q=add(k,neg(p))
            if q not in signed:continue
            fac=wk*cross(p,q)*Q(1,norm2(p))
            if not fac:continue
            row=phase_row(modes,[(p,1),(q,1),(k,-1)])
            fac*=amps[canonical(p)]*amps[canonical(q)]*amps[canonical(k)]
            if row:out[row]=out.get(row,Q(0))+fac
    return {r:out[r] for r in sorted(out) if out[r]}

def evaluate(coeff, phases):
    return sum(float(c)*np.cos(np.dot(r,phases)) for r,c in coeff.items())

def four_cycle_max(weights):
    """Sharp max sum c_i cos(theta_i), theta1+theta3=theta2+theta4.

    Returns rational maximizer t=cos(theta1+theta3), squared radicands,
    and a floating display. The rational quantities specify the exact answer.
    """
    c1,c2,c3,c4=map(Q,weights)
    a=c1*c1+c3*c3;b=2*c1*c3;c=c2*c2+c4*c4;d=2*c2*c4
    candidates=[Q(-1),Q(1)]
    if b*d<0:
        t=(d*d*a-b*b*c)/(b*b*d-d*d*b)
        if -1<=t<=1:candidates.append(t)
    def f(t):return np.sqrt(float(a+b*t))+np.sqrt(float(c+d*t))
    # With bd<0 the derivative is strictly decreasing; an interior stationary
    # point, when present, is the global maximum. Otherwise monotonicity or
    # concavity puts the maximum at an endpoint. Endpoint comparisons are exact.
    if len(candidates)==3:t=candidates[-1]
    else:
        # sqrt(a +/- b)=|c1 +/- c3|, similarly for second pair.
        fm=abs(c1-c3)+abs(c2-c4);fp=abs(c1+c3)+abs(c2+c4)
        t=Q(1) if fp>=fm else Q(-1)
    return {'t':t,'radicands':(a+b*t,c+d*t),'maximum':f(t),
            'triangle':sum(map(abs,(c1,c2,c3,c4)))}

def maximizing_phases(weights,t):
    """Attaining cycle phases (display precision), constructively by phasors."""
    c1,c2,c3,c4=map(float,weights);s=np.arccos(float(t))
    # theta3=s-theta1; theta4=s-theta2.
    x=-np.angle(c1+c3*np.exp(-1j*s))
    z=-np.angle(c2+c4*np.exp(-1j*s))
    return np.array([x,z,s-x,s-z])

def four_circuits(rows):
    """Enumerate distinct four-row relations with coefficients +/-1.

    Pair-sum hash, including signed rows; canonical relation starts positive.
    This finds four-circuits only, not every integer dependence.
    """
    rows=list(rows);pairs={};found=set()
    for i,j in combinations(range(len(rows)),2):
        for sign in (-1,1):
            v=tuple(a+sign*b for a,b in zip(rows[i],rows[j]))
            orient=1 if next((x for x in v if x),1)>0 else -1
            key=tuple(orient*x for x in v)
            entry=(i,j,orient,orient*sign)
            for k,l,sk,sl in pairs.get(key,[]):
                if len({i,j,k,l})!=4:continue
                rel=sorted([(i,orient),(j,orient*sign),(k,-sk),(l,-sl)])
                if rel[0][1]<0:rel=[(a,-b) for a,b in rel]
                found.add(tuple(rel))
            pairs.setdefault(key,[]).append(entry)
    return sorted(found)

def pack_certificate(coeff):
    """Rational feasible equal-weight circuit packing; no optimizer trust.

    Valid upper bound U=triangle-(4-2sqrt(2))*packed_mass.
    Greedy order uses descending initial minimum weight then lexicographic order.
    """
    rows=list(coeff);cs=list(coeff.values());remaining=list(map(abs,cs));cycles=[]
    for rel in four_circuits(rows):
        ids=[i for i,_ in rel]
        if np.prod([1 if cs[i]>0 else -1 for i in ids])==-1:
            cycles.append(rel)
    cycles.sort(key=lambda rel:(-min(abs(cs[i]) for i,_ in rel),rel))
    packed=[]
    for rel in cycles:
        mass=min(remaining[i] for i,_ in rel)
        if mass:
            for i,_ in rel:remaining[i]-=mass
            packed.append((rel,mass))
    mass=sum((m for _,m in packed),Q(0));tri=sum(map(abs,cs),Q(0))
    return {'triangle':tri,'packed_mass':mass,'cycles':packed,
            'upper':float(tri)-(4-2*np.sqrt(2))*float(mass),
            'candidate_cycles':len(cycles)}

def verify_packing(coeff, cert):
    rows=list(coeff);cs=list(coeff.values());used=[Q(0)]*len(rows)
    for rel,m in cert['cycles']:
        assert m>0 and len(set(i for i,_ in rel))==4
        assert all(s in (-1,1) for _,s in rel)
        assert all(sum(s*rows[i][j] for i,s in rel)==0 for j in range(len(rows[0])))
        assert np.prod([1 if cs[i]>0 else -1 for i,_ in rel])==-1
        for i,_ in rel:used[i]+=m
    assert all(u<=abs(c) for u,c in zip(used,cs))
    assert sum((m for _,m in cert['cycles']),Q(0))==cert['packed_mass']
    assert sum(map(abs,cs),Q(0))==cert['triangle']
    return True

WITNESS_MODES=[(1,-3),(1,-2),(1,-1),(2,-5),(2,-3),(3,-6)]
WITNESS_ROWS=[(1,1,0,-1,0,0),(0,1,1,0,-1,0),(0,0,1,1,0,-1),(1,0,0,0,1,-1)]

def integer_relations(rows):
    """Primitive exact integer relations from a rational left-nullspace basis.

    These are sufficient witnesses; no claim that the returned primitive vectors
    form a saturated integer-lattice basis is needed for certificate validity.
    """
    import sympy as sp
    if not rows:return []
    out=[]
    for v in sp.Matrix(list(rows)).T.nullspace():
        denominator=1
        for x in v:denominator=sp.ilcm(denominator,x.q)
        z=[int(x*denominator) for x in v];g=0
        for x in z:g=gcd(g,abs(x))
        z=[x//g for x in z]
        if next(x for x in z if x)<0:z=[-x for x in z]
        out.append(tuple((i,x) for i,x in enumerate(z) if x))
    return out

def pack_hybrid_certificate(coeff):
    """Four-row sharp penalties followed by exact general-relation penalties.

    For an integer relation z and equally allocated mass m on its support,
    sign incompatibility forces deficit >= 2*m/sum(z_i^2).
    """
    base=pack_certificate(coeff);rows=list(coeff);cs=list(coeff.values())
    remaining=list(map(abs,cs))
    for rel,m in base['cycles']:
        for i,_ in rel:remaining[i]-=m
    relations=[]
    for rel in integer_relations(rows):
        parity=sum(z for i,z in rel if cs[i]<0)%2
        if parity:
            sq=sum(z*z for _,z in rel)
            relations.append((rel,Q(2,sq)))
    relations.sort(key=lambda x:(-x[1]*min(remaining[i] for i,_ in x[0]),x[0]))
    packed=[]
    for rel,penalty in relations:
        mass=min(remaining[i] for i,_ in rel)
        if mass:
            for i,_ in rel:remaining[i]-=mass
            packed.append((rel,mass,penalty))
    deficit=sum((m*p for _,m,p in packed),Q(0))
    return {**base,'general_relations':packed,'general_deficit':deficit,
            'general_candidates':len(relations),'upper':base['upper']-float(deficit)}

def verify_hybrid(coeff,cert):
    verify_packing(coeff,cert);rows=list(coeff);cs=list(coeff.values());used=[Q(0)]*len(rows)
    for rel,m in cert['cycles']:
        for i,_ in rel:used[i]+=m
    deficit=Q(0)
    for rel,m,penalty in cert['general_relations']:
        assert m>0
        assert all(sum(z*rows[i][j] for i,z in rel)==0 for j in range(len(rows[0])))
        assert sum(z for i,z in rel if cs[i]<0)%2==1
        assert penalty==Q(2,sum(z*z for _,z in rel))
        deficit+=m*penalty
        for i,_ in rel:used[i]+=m
    assert all(u<=abs(c) for u,c in zip(used,cs))
    assert deficit==cert['general_deficit']
    return True

def triangle_attainable(rows,coeff):
    """Complete attainability test using integer Smith decomposition.

    Multiplication by a nonzero integer is surjective on the continuous phase
    torus; testing Mx=b over F2 instead would incorrectly reject torsion cases.
    """
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix
    from sympy.polys.matrices.normalforms import smith_normal_decomp
    if not rows:return True
    if any(not c for c in coeff):raise ValueError('remove zero coefficients first')
    M=DomainMatrix.from_Matrix(sp.Matrix(rows)).convert_to(sp.ZZ)
    D,U,V=smith_normal_decomp(M)
    dd=D.to_Matrix();u=U.to_Matrix();rank=sum(dd[i,i]!=0 for i in range(min(dd.shape)))
    b=sp.Matrix([int(c<0) for c in coeff]);target=u*b
    return all(int(target[i])%2==0 for i in range(rank,len(coeff)))
