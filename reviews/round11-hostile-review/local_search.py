import math,json,copy,itertools
from pathlib import Path
from scipy.optimize import differential_evolution
from review_checks import load_module
A=load_module('upload/appd_assembly.py')
C=json.load(open('upload/certificates.json'))['certificates']

def build(rs,us,delta,Lc,kap=10):
 p=dict(levels=len(rs),r=rs,u=us,delta=delta,Lc=Lc,a=[],kappa=[kap]*len(rs),
        **{f:True for f in A.FLAGS})
 p['leafvars']=len(rs)==1
 M=1+1/Lc;lp=0.
 for r,u in zip(rs,us):
  aa=max(u+kap,math.ceil(lp-math.log2(M)+u+1))
  while A.need_bits(lp,M,aa,u,kap)>aa-u+1:aa+=1
  while aa>u+kap and A.need_bits(lp,M,aa-1,u,kap)<=aa-u:aa-=1
  p['a'].append(aa)
  lp=A.next_lp(lp,M,aa,u,kap)
  M*=2*(1+delta)/r
 return p

def required(p,matching=False):
 M=1+1/p['Lc']; req=64.;by=[]
 for d,(r,u,a) in enumerate(zip(p['r'],p['u'],p['a'])):
  D=2*u*math.log2(r); lM=math.log2(M)
  if d==0:
   if matching:
    q=126*p['Lc']-1
    slack=math.log2(8*q)+2*math.log2(1/p['delta'])+math.log2(math.log(8*q)+D*math.log(2))
   else:
    slack=13+2*math.log2(p['Lc'])+2*math.log2(1/p['delta'])+math.log2(math.log(172*p['Lc'])+D*math.log(2))
  else:slack=3+p['a'][d-1]+1+2*math.log2(1/p['delta'])+math.log2(D*math.log(2))
  bal=D+slack-lM
  pad=a+u*math.log2(r*(r+1)/2)-lM-math.log2(p['delta'])
  req=max(req,bal,pad);by.append([bal,pad])
  M*=2*(1+p['delta'])/r
 return req,by

out={}
for matching in (False,True):
 for name,c in C.items():
  p0=c['params']; L=p0['levels'];tar=A.TARGETS[name]
  best=[float('inf'),None,None]
  def fun(x):
   rs=[int(round(y)) for y in x[:L]]
   us=[int(round(y)) for y in x[L:2*L]]
   delta=10**x[-2];Lc=2**int(round(x[-1]))
   p=build(rs,us,delta,Lc)
   req,by=required(p,matching)
   # Frozen code at a deliberately generous lambda, so only the cost is queried.
   val=A.exponent_of(p,100000)
   if not val['valid']:return 1e6
   e=val['e']
   if e<tar and req<best[0]:best[:]=[req,p,e]
   return req+max(e/tar-1,0)*40000
  bounds=[(max(4,r-3),min(16,r+3)) for r in p0['r']]
  bounds +=[(max(2,int(u*.65)),int(u*1.4)+2) for u in p0['u']]
  bounds +=[(-4,-.65),(3,20)]
  x0=p0['r']+p0['u']+[math.log10(p0['delta']),math.log2(p0['Lc'])]
  fun(x0)
  res=differential_evolution(fun,bounds,seed=72,maxiter=160,popsize=12,polish=False,x0=x0,tol=1e-8)
  req,p,e=best
  key=name+('_matching' if matching else '')
  out[key]=dict(required_lambda=req,lam=math.ceil(req),params=p,e=e,checks=required(p,matching)[1])
  if not matching:
   check=A.exponent_of(p,math.ceil(req));assert check['valid'] and check['e']<tar
  print(key,round(req,6),math.ceil(req),e,flush=True)
  Path('search_results.json').write_text(json.dumps(out,indent=2)+'\n')
