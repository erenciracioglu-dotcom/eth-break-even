"""Independent arithmetic checks for the supplied Appendix D model.

Run: python3 review_checks.py --assembly upload/appd_assembly.py \
       --certificates upload/certificates.json --output review_checks.json
Requires numpy for batched exact-integer profile certificates.
This is not a reimplementation of the omitted SAT counter or referee.
"""
import argparse
import importlib.util
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path


def load_module(path):
    spec = importlib.util.spec_from_file_location('assembly', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def independent_forms(order):
    r = len(order)
    edges = list(itertools.combinations(range(r), 2))
    forms = set()
    for split in range(1, r):
        A, B = set(order[:split]), {order[split]}
        C = set(order[split+1:])
        cuts = [(A, B), (A, C), (B, C)]
        for discounted in range(3):
            f = []
            for v, w in edges:
                val = 0
                for k, (left, right) in enumerate(cuts):
                    if (v in left and w in right) or (w in left and v in right):
                        val = 101 if k == discounted else 125
                f.append(val)
            forms.add(tuple(f))
    return sorted(forms)


def independent_certificate(menu):
    import numpy as np
    r, count = len(menu[0]), len(menu)
    arrays = [np.array(independent_forms(o), dtype=np.int64) for o in menu]
    # At most four forms, each coordinate <= 125. Sum of squares < 2^31.
    left = arrays[0]
    for ar in arrays[1:-1]:
        left = (left[:, None, :] + ar[None, :, :]).reshape(-1, left.shape[1])
    best = 0
    for start in range(0, len(left), 256):
        sums = left[start:start+256, None, :] + arrays[-1][None, :, :]
        best = max(best, int(np.max(np.sum(sums*sums, axis=-1))))
    return F(2*best, count*count*125*125*r*r)


def v2(x):
    assert x
    return (abs(x) & -abs(x)).bit_length() - 1


def check_hermite_small():
    # Exact valuation of the cardinal factor, directly from its product.
    cases = 0
    for b in range(2, 9):
        H = 1 << b
        for u in range(2, b+1):
            q = {t: max(0, u-v2(t)) for t in range(-H, H) if t}
            for t in q:
                i = v2(t)
                if not q[t]:
                    continue
                direct = sum(mult*(v2(t-s)-v2(s))
                             for s, mult in q.items() if s != t)
                formula = sum(F(H*j*min(j, u-i), 1 << (i+j))
                              for j in range(1, b-i+1))
                assert direct == formula
                assert 6*F(H, 1 << i) - direct >= 4*F(H, 1 << (u-1))
                cases += 1
    return cases


def finite_majorization():
    cases = 0
    for k in range(17):
        ws = [k+2+4*b-2*b.bit_count() for b in range(1 << k)]
        assert ws == sorted(ws)
        pw = pv = 0
        for b, w in enumerate(ws):
            pw += w
            pv += 4*b+2
            assert pw >= pv
        assert pw == pv
        cases += 1
    return cases


def binary_factorizations():
    import numpy as np
    I=np.eye(4,dtype=np.int64)
    v=np.array([0,1,-1,0],dtype=np.int64)
    J=np.outer(v,v)
    swap=I-J
    assert np.array_equal(swap,np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]))
    def colflat(x):
        return x.reshape(2,2,2,2).transpose(0,2,1,3).reshape(4,4)
    K=colflat(J);P0=colflat(I)
    assert np.array_equal(K@K,I) and np.array_equal(K@P0,P0)
    S=np.array([[0,1,0,0],[1,0,0,0],[-1,0,1,0],[0,0,0,1]],dtype=np.int64)
    T=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[1,0,0,1]],dtype=np.int64)
    Sinv=np.rint(np.linalg.inv(S)).astype(np.int64)
    Tinv=np.rint(np.linalg.inv(T)).astype(np.int64)
    assert np.array_equal(S@Sinv,I) and np.array_equal(T@Tinv,I)
    def power(x,m):
        ans=np.ones((1,1),dtype=np.int64)
        for _ in range(m):ans=np.kron(ans,x)
        return ans
    def coeffs(A,B,m):
        c=[np.ones((1,1),dtype=np.int64)]
        for _ in range(m):
            new=[np.zeros((4*c[0].shape[0],)*2,dtype=np.int64) for j in range(len(c)+1)]
            for j,x in enumerate(c):
                new[j]+=np.kron(x,A);new[j+1]+=np.kron(x,B)
            c=new
        return c
    cases=0
    for m in range(5):
        zeros=[sum(d==0 for d in word) for word in itertools.product(range(4),repeat=m)]
        A=coeffs(I,-J,m);B=coeffs(P0,-I,m)
        assert np.array_equal(sum(A),power(swap,m))
        for i in range(m+1):
            R=power(Sinv,m)@A[i]
            mask=np.array([z>=i for z in zeros])
            assert not np.any(R[~mask])
            assert np.array_equal(power(S,m)[:,mask]@R[mask],A[i])
            C=power(Tinv,m)@B[i]
            mask=np.array([z>=m-i for z in zeros])
            assert not np.any(C[~mask])
            assert np.array_equal(power(K@T,m)[:,mask]@C[mask],power(K,m)@B[i])
            cases+=1
    return cases


def precision_certificate(params, lam):
    # All p,M,h and losses below are exact rationals/integers. Test every
    # retained valuation i, using the complete residue-grid formula.
    n = 1 << lam
    delta = F(str(params['delta']))
    M = n*(1+F(1, params['Lc']))
    pin = n+1
    results = []
    for r,u,a,kap in zip(params['r'],params['u'],params['a'],params['kappa']):
        Q,H = 1 << a, 1 << (a-kap)
        # Take the integer wire count up to its stated rational bound.
        Mi = (M.numerator+M.denominator-1)//M.denominator
        S = sum((H >> i)*(u-i) for i in range(u))
        degree = Mi*(Q-1)
        h = degree//S+1
        # Frozen model's K*p uses the following effective precision.
        # h_model * 2(Q-1)/(2Q-1), including the model's guarded ratio.
        h_eff = F(2**(kap-1)*Mi, u-1)*F(2*(Q-1), 2*Q-1)
        min_relative = None
        for i in range(u):
            D = F(H, 1 << i)
            loss = h*sum(F(H*j*min(j,u-i), 1 << (i+j))
                         for j in range(1, a-kap-i+1))
            required = pin+loss+(h*(u-i)-1)*(a+1-i)
            allowed = 6*h_eff*D
            rel = (allowed-required)/allowed
            min_relative = rel if min_relative is None else min(min_relative,rel)
            assert required <= allowed, (lam,r,u,a,i)
        assert h*S > degree
        results.append({'r':r,'u':u,'a':a,'min_relative_precision_slack':float(min_relative)})
        # Conservative inherited precision charged by next_lp, in exact arithmetic.
        pin = pin+F(3*Mi*(1 << a),u-1)+(1 << (kap+1))*Mi*(a+1)
        M *= F(2)*(1+delta)/r
    return results


def matching_size_conditions(params, lam):
    """Proposed replacement only for the root balance check; see review proof.

    q <= 2*Delta-1, Delta <= 63*Lc. The sufficient load condition is
    M/D >= 8*q*log(8*q*D)/delta^2. All other frozen checks are retained.
    """
    M=1+1/params['Lc']
    requirements=[]
    for depth,(r,u,a) in enumerate(zip(params['r'],params['u'],params['a'])):
        logD=2*u*math.log2(r)
        if depth==0:
            q=126*params['Lc']-1
            slack=math.log2(8*q)+2*math.log2(1/params['delta'])+math.log2(math.log(8*q)+logD*math.log(2))
        else:
            slack=4+params['a'][depth-1]+2*math.log2(1/params['delta'])+math.log2(logD*math.log(2))
        balance=logD+slack-math.log2(M)
        padding=a+u*math.log2(r*(r+1)/2)-math.log2(M)-math.log2(params['delta'])
        requirements.append({'depth':depth,'balance_required_lambda':balance,
                             'padding_required_lambda':padding})
        assert max(balance,padding) <= lam
        M *= 2*(1+params['delta'])/r
    return requirements


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--assembly',default='upload/appd_assembly.py')
    ap.add_argument('--certificates',default='upload/certificates.json')
    ap.add_argument('--output',default='review_checks.json')
    ap.add_argument('--additional',default='additional_certificates.json')
    args=ap.parse_args()
    A=load_module(args.assembly)
    supplied=json.loads(Path(args.certificates).read_text())
    out={'scope':'Arithmetic checks; no original algorithm/referee source supplied',
         'certificates':{},'profile_certificates':{}}
    for name,c in supplied['certificates'].items():
        got=A.exponent_of(c['params'],c['lam'])
        assert got['valid'] and got['e']==c['e']
        got['margin']=A.TARGETS[name]-got['e']
        got['previous_size']=A.exponent_of(c['params'],c['lam']-1)
        got['exact_precision_checks']=precision_certificate(c['params'],c['lam'])
        out['certificates'][name]=got
    for r,menu in A.PF_MENUS.items():
        c=independent_certificate(menu)
        assert c==A.PF_EPS[r], (r,c,A.PF_EPS[r])
        out['profile_certificates'][r]=str(c)
    assert 7**125 < 2**351
    out['hermite_grid_cases']=check_hermite_small()
    out['finite_ldu_majorization_cases']=finite_majorization()
    out['binary_factorization_pieces']=binary_factorizations()
    out['binary_channel_cases']=0
    for m in range(101):
        bound=sum(math.comb(m,j)*3**j for j in range(m//2+1))
        assert bound*bound <= 12**m
        out['binary_channel_cases']+=1
    if Path(args.additional).exists():
        out['additional_certificates']={}
        extra=json.loads(Path(args.additional).read_text())
        for key,c in extra.items():
            target=key.split('_')[0]
            if c['model']=='frozen-unchanged':
                got=A.exponent_of(c['params'],c['lam'])
            else:
                # Cost only from frozen model; size feasibility from the
                # proposed matching lemma, never claim frozen acceptance.
                got=A.exponent_of(c['params'],100000)
                c['matching_size_checks']=matching_size_conditions(c['params'],c['lam'])
            assert got['valid'] and got['e'] < A.TARGETS[target]
            assert got['e']==c['e']
            c['exact_precision_checks']=precision_certificate(c['params'],c['lam'])
            out['additional_certificates'][key]=c
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'profiles':len(out['profile_certificates']),
                      'hermite_cases':out['hermite_grid_cases'],
                      'output':args.output}))


if __name__=='__main__':
    main()
