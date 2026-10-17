"""Reproduce all paper figures, tables, exact certificate data and numerical checks."""
import sys,json,time,csv,platform
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from triad_flux import *
from independent import state,flux,direct_flux,rk4
import scipy
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sympy as sp
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';FIG=ROOT/'paper/figures'
OUT.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
rng=np.random.default_rng(20261017)
def save(name):
    plt.savefig(FIG/(name+'.pdf'));plt.savefig(FIG/(name+'.png'),dpi=180);plt.close()
def dump(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,default=lambda x:str(x))+'\n')

# Exact witness, sharp range, explicit attaining modes.
C=coefficients(WITNESS_MODES,10);weights=[C[r] for r in WITNESS_ROWS]
sharp=four_cycle_max(weights);theta=maximizing_phases(weights,sharp['t'])
phi=np.linalg.lstsq(np.array(WITNESS_ROWS,float),theta,rcond=None)[0]
assert np.max(abs(np.array(WITNESS_ROWS)@phi-theta))<1e-12
assert abs(evaluate(C,phi)-sharp['maximum'])<1e-13
# Algebraic balanced amplitudes, with weights (+,+,-,+)/145 exactly.
k1,k2,k3,k4=map(sp.Rational,weights)
c=sp.sqrt(k4*k1/(k2*(-k3)))
amp=[sp.Integer(1),sp.Integer(1),c,sp.Integer(1),sp.simplify(k1/(k2*c)),sp.simplify(k1/((-k3)*c))]
products=[amp[0]*amp[1]*amp[3],amp[1]*amp[2]*amp[4],amp[3]*amp[2]*amp[5],amp[0]*amp[4]*amp[5]]
assert [sp.simplify(k*p) for k,p in zip(map(sp.Rational,weights),products)]==[k1,k1,-k1,k1]
balanced=four_cycle_max([Q(1,145),Q(1,145),-Q(1,145),Q(1,145)])
btheta=maximizing_phases([1,1,-1,1],balanced['t'])
bphi=np.linalg.lstsq(np.array(WITNESS_ROWS,float),btheta,rcond=None)[0]
# Independent representations, 600 evaluations in 3 grids and 2 observables.
errors=[]
for n in [40,48,64]:
 for obs in ['energy','enstrophy']:
  for seed in range(100):
   phases=rng.uniform(-np.pi,np.pi,6);amps=[Q(int(i),13) for i in rng.integers(1,25,6)]
   coeff=coefficients(WITNESS_MODES,10,obs,amps)
   grouped=evaluate(coeff,phases);direct=direct_flux(WITNESS_MODES,amps,phases,10,obs)
   fft=flux(state(WITNESS_MODES,np.array(amps,float),phases,n),10,obs)
   errors.append([n,obs,seed,abs(grouped-direct),abs(grouped-fft)])
with (OUT/'representation_checks.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['grid','observable','replicate','direct_error','fft_error']);w.writerows(errors)

# Complete declared weighted-cycle grid, 5^4 - zero vector = 624 cases.
# Each exact value accompanied by an explicit phase witness; optimizer is a diagnostic.
cycle_tests=[]
from itertools import product
for ctuple in product([-2,-1,0,1,2],repeat=4):
 if not any(ctuple):continue
 ans=four_cycle_max(ctuple);th=maximizing_phases(ctuple,ans['t'])
 attained=float(np.array(ctuple)@np.cos(th))
 assert abs(attained-ans['maximum'])<2e-12
 # Concave reduction compared to a dense scan, not used as certification.
 tt=np.linspace(-1,1,1001);c1,c2,c3,c4=ctuple
 scan=np.sqrt(np.maximum(0,c1*c1+c3*c3+2*c1*c3*tt))+np.sqrt(np.maximum(0,c2*c2+c4*c4+2*c2*c4*tt))
 assert max(scan)<=ans['maximum']+2e-12
 cycle_tests.append([*ctuple,str(ans['t']),ans['maximum'],abs(attained-ans['maximum']),ans['maximum']-max(scan)])
with (OUT/'cycle_checks.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['c1','c2','c3','c4','t_exact','sharp_max','attainment_error','scan_shortfall']);w.writerows(cycle_tests)

# Dense finite supports: every mode in [-L,L]^2, exact rational amplitudes.
# All declared conditions retained, including zero-improvement conditions.
dense=[];certificates=[]
for L in [2,3,4]:
 modes=[(x,y) for x in range(L+1) for y in range(-L,L+1) if x>0 or (x==0 and y>0)]
 for decay in [0,1]:
  amps=[Q(1,(1+norm2(k))**decay) for k in modes]
  for obs in ['energy','enstrophy']:
   cutoff2=L*L
   start=time.perf_counter();cc=coefficients(modes,cutoff2,obs,amps);cert=pack_hybrid_certificate(cc);elapsed=time.perf_counter()-start
   verify_hybrid(cc,cert)
   rows=np.array(list(cc),float);cs=np.array(list(cc.values()),float)
   fun=lambda p:-np.dot(cs,np.cos(rows@p))
   jac=lambda p: rows.T@(cs*np.sin(rows@p))
   lower=-np.inf
   for seed in range(12):
    p=rng.uniform(-np.pi,np.pi,len(modes))
    opt=minimize(fun,p,jac=jac,method='L-BFGS-B',options={'maxiter':1200,'ftol':1e-13,'gtol':1e-9})
    lower=max(lower,-opt.fun)
   assert lower<=cert['upper']+1e-9*max(1,float(cert['triangle']))
   record={'L':L,'modes':len(modes),'decay':decay,'observable':obs,'triads':len(cc),'triangle_attainable':bool(triangle_attainable(list(cc),list(cc.values()))),'frustrated_circuits':cert['candidate_cycles'],'used_circuits':len(cert['cycles']),'general_candidates':cert['general_candidates'],'used_general':len(cert['general_relations']),'general_deficit':float(cert['general_deficit']),'triangle':float(cert['triangle']),'upper':cert['upper'],'lower_diagnostic':lower,'improvement_percent':100*(1-cert['upper']/float(cert['triangle'])),'certificate_seconds':elapsed}
   dense.append(record)
   certificates.append({'condition':record,'modes':modes,'amplitudes':list(map(str,amps)),'cutoff2':cutoff2,'rows':list(cc),'coefficients':list(map(str,cc.values())),'triangle':str(cert['triangle']),'packed_mass':str(cert['packed_mass']),'cycles':[{'relation':r,'mass':str(m)} for r,m in cert['cycles']],'general_relations':[{'relation':r,'mass':str(m),'penalty':str(p)} for r,m,p in cert['general_relations']],'general_deficit':str(cert['general_deficit'])})
   print(record,flush=True)
dump('packing_certificates.json',certificates)
with (OUT/'dense_results.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(dense[0]));w.writeheader();w.writerows(dense)

# Short-time NS trajectories from the balanced maximizing initial state.
# These are resolution checks only; instantaneous certificate is analytic.
trajectories=[];nu=.01;T=.2
for n,dt in [(48,.002),(64,.002),(64,.001),(96,.001)]:
 w=state(WITNESS_MODES,np.array(amp,float),bphi,n)
 initial=flux(w,10);track=[]
 k=np.rint(np.fft.fftfreq(n)*n);r=k[:,None]**2+k[None,:]**2
 for step in range(round(T/dt)+1):
  if step%round(.01/dt)==0:
   E=.5*np.sum(abs(w)**2/np.where(r==0,1,r));Z=.5*np.sum(abs(w)**2)
   track.append([step*dt,flux(w,10),E,Z])
  if step<round(T/dt):w=rk4(w,dt,nu)
 trajectories.append({'n':n,'dt':dt,'viscosity':nu,'track':track})
dump('trajectories.json',trajectories)
# Forward finite difference checks of shell-energy derivative (nonlinear plus viscosity).
diffs=[]
for dt in [.01,.005,.0025,.00125,.000625]:
 w=state(WITNESS_MODES,np.array(amp,float),bphi,64);k=np.rint(np.fft.fftfreq(64)*64);r=k[:,None]**2+k[None,:]**2;mask=(r>0)&(r<=10)
 shell=lambda v:.5*np.sum(abs(v[mask])**2/r[mask])
 loss=-(shell(rk4(w,dt,nu))-shell(w))/dt
 target=balanced['maximum']+nu*np.sum(abs(w[mask])**2)
 diffs.append([dt,loss,target,abs(loss-target)])
dump('derivative_checks.json',diffs)

# Figure 1: Fourier support and cycle relation, with physical vorticity witness.
fig,ax=plt.subplots(1,2,figsize=(7,3.1))
S=np.array(WITNESS_MODES)
ax[0].scatter(S[:,0],S[:,1],s=60,c=['#007c91']*3+['#ce6b24']*3)
ax[0].scatter(-S[:,0],-S[:,1],s=30,color='#aab7bd')
for i,p in enumerate(S):ax[0].annotate('abcdef'[i],p,xytext=(5,4),textcoords='offset points')
ang=np.linspace(0,2*np.pi,300);ax[0].plot(np.sqrt(10)*np.cos(ang),np.sqrt(10)*np.sin(ang),'--',color='.5',lw=1)
ax[0].set(xlabel='$k_x$',ylabel='$k_y$',title='(a) Exact occupied Fourier modes');ax[0].set_aspect('equal')
w=state(WITNESS_MODES,np.array(amp,float),bphi,128);v=np.fft.ifft2(w*128**2).real
im=ax[1].imshow(v.T,origin='lower',extent=[0,2*np.pi,0,2*np.pi],cmap='RdBu_r');ax[1].set(xlabel='$x$',ylabel='$y$',title='(b) Attaining vorticity field');fig.colorbar(im,ax=ax[1],shrink=.75,label='$\omega$')
fig.tight_layout();save('support_witness')

# Figure 2: exact 1D objective for unit amplitudes and balanced amplitudes.
fig,ax=plt.subplots(1,2,figsize=(7,2.8));tt=np.linspace(-1,1,1001)
for axis,ww,title in zip(ax,[weights,[1,1,-1,1]],['Unit vorticity amplitudes','Balanced transfer weights']):
 c1,c2,c3,c4=map(float,ww);tri=sum(abs(np.array(ww,float)))
 ff=np.sqrt(np.maximum(0,c1*c1+c3*c3+2*c1*c3*tt))+np.sqrt(np.maximum(0,c2*c2+c4*c4+2*c2*c4*tt))
 axis.plot(tt,ff/tri,color='#007c91',lw=2);axis.axhline(1,color='#b63b36',ls='--',label='Independent-triad bound')
 a=four_cycle_max(ww);axis.scatter([float(a['t'])],[a['maximum']/tri],c='#ce6b24',zorder=3)
 axis.set(xlabel=r'$t=\cos(\theta_1+\theta_3)$',ylabel='Flux / triangle bound',title=title,ylim=(0,1.07));axis.grid(alpha=.2)
ax[1].legend(fontsize=8,loc='lower center');fig.tight_layout();save('sharp_envelopes')

# Figure 3: dense-support bound tightening, and unattested optimizer lower values.
fig,ax=plt.subplots(1,2,figsize=(7,3))
for obs,axis in zip(['energy','enstrophy'],ax):
 for decay,marker in [(0,'o'),(1,'s')]:
  rr=[r for r in dense if r['observable']==obs and r['decay']==decay]
  axis.plot([r['L'] for r in rr],[r['improvement_percent'] for r in rr],marker=marker,label=f'amplitude $(1+|k|^2)^{{-{decay}}}$')
 axis.set(xlabel='Support half-width $L$',ylabel='Reduction of upper bound (%)',title=obs.capitalize());axis.set_xticks([2,3,4]);axis.grid(alpha=.2)
ax[1].legend(fontsize=8);fig.tight_layout();save('packing_gains')

# Figure 4: short-time resolution and derivative consistency.
fig,ax=plt.subplots(1,2,figsize=(7,2.9))
for tr in trajectories:
 t,f,*_=np.array(tr['track']).T;ax[0].plot(t,f,label=f"$N={tr['n']}$, $h={tr['dt']}$",lw=1.4)
ax[0].scatter([0],[balanced['maximum']],color='black',s=16,zorder=3,label='Exact initial maximum')
ax[0].set(xlabel='Time',ylabel='Nonlinear shell-energy flux',title=r'(a) Navier--Stokes, $\nu=0.01$');ax[0].legend(fontsize=7)
d=np.array(diffs);ax[1].loglog(d[:,0],d[:,3],'o-',color='#007c91');ax[1].loglog(d[:,0],d[-1,3]*d[:,0]/d[-1,0],'--',color='.5',label='First-order reference');ax[1].set(xlabel='Forward difference interval',ylabel='Derivative error',title='(b) Shell-energy consistency');ax[1].legend(fontsize=8);fig.tight_layout();save('time_checks')

summary={'witness_weights':list(map(str,weights)),'witness_sharp':sharp,'witness_phases':phi.tolist(),'balanced_amplitudes':list(map(str,amp)),'balanced_phases':bphi.tolist(),'balanced_sharp':balanced,'representation_count':len(errors),'max_direct_error':max(x[3] for x in errors),'max_fft_error':max(x[4] for x in errors),'cycle_count':len(cycle_tests),'max_cycle_attainment_error':max(x[6] for x in cycle_tests),'trajectory_terminal_fluxes':[x['track'][-1][1] for x in trajectories],'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,'matplotlib':matplotlib.__version__,'seed':20261017}
dump('summary.json',summary)
print(json.dumps(summary,indent=2,default=str))
