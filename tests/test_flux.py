import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from triad_flux import *
from independent import state,flux,direct_flux,rhs
from scipy.optimize import differential_evolution

class FluxTests(unittest.TestCase):
    def test_witness_coefficients(self):
        c=coefficients(WITNESS_MODES,10)
        self.assertEqual([c[r] for r in WITNESS_ROWS],[Q(1,145),Q(3,65),-Q(16,435),Q(32,975)])
        self.assertEqual(len(c),4)
        self.assertEqual(tuple(np.array(WITNESS_ROWS).T@np.array([1,-1,1,-1])),(0,)*6)
    def test_full_invariant_cancellation(self):
        for L in (1,2,3):
            modes=[(x,y) for x in range(L+1) for y in range(-L,L+1) if x>0 or (x==0 and y>0)]
            for obs in ('energy','enstrophy'):self.assertEqual(coefficients(modes,2*L*L,obs),{})
    def test_three_evaluators(self):
        rng=np.random.default_rng(63419)
        for n in (40,48,64):
            for obs in ('energy','enstrophy'):
                for _ in range(25):
                    phi=rng.uniform(-np.pi,np.pi,6);a=[Q(int(x),7) for x in rng.integers(1,15,6)]
                    c=coefficients(WITNESS_MODES,10,obs,a)
                    f=evaluate(c,phi);g=direct_flux(WITNESS_MODES,a,phi,10,obs)
                    h=flux(state(WITNESS_MODES,np.array(a,float),phi,n),10,obs)
                    self.assertAlmostEqual(f,g,places=11);self.assertAlmostEqual(f,h,places=11)
    def test_cycle_attainment(self):
        rng=np.random.default_rng(59310)
        for _ in range(160):
            c=[Q(int(x),17) for x in rng.integers(-20,21,4)]
            ans=four_cycle_max(c);theta=maximizing_phases(c,ans['t'])
            self.assertAlmostEqual(theta[0]+theta[2],theta[1]+theta[3],places=12)
            f=sum(float(x)*np.cos(y) for x,y in zip(c,theta))
            self.assertAlmostEqual(f,ans['maximum'],places=11)
            random_theta=rng.uniform(-np.pi,np.pi,(200,3));th=np.column_stack([random_theta,random_theta[:,0]+random_theta[:,2]-random_theta[:,1]])
            self.assertLessEqual(np.max(np.cos(th)@np.array(c,float)),ans['maximum']+1e-12)
    def test_balanced_frustration(self):
        for i in range(16):
            c=[Q((-1)**((i>>j)&1)) for j in range(4)]
            ans=four_cycle_max(c)
            expected=2*np.sqrt(2) if np.prod(c)==-1 else 4
            self.assertAlmostEqual(ans['maximum'],expected,places=12)
    def test_packing_and_reject_tamper(self):
        modes=[(x,y) for x in range(3) for y in range(-2,3) if x>0 or (x==0 and y>0)]
        c=coefficients(WITNESS_MODES,10);cert=pack_certificate(c)
        self.assertTrue(verify_packing(c,cert))
        self.assertGreater(cert['packed_mass'],0)
        old=cert['packed_mass'];cert['packed_mass']=old+1
        with self.assertRaises(AssertionError):verify_packing(c,cert)
    def test_viscous_split(self):
        phi=np.arange(6)*.21;w=state(WITNESS_MODES,np.ones(6),phi,48)
        k=np.fft.fftfreq(48)*48;r=k[:,None]**2+k[None,:]**2;mask=(r>0)&(r<=10)
        nu=.03;weight=1/np.where(r==0,1,r)
        total=-np.sum((weight*w.conj()*rhs(w,nu)).real[mask])
        decay=nu*np.sum(abs(w[mask])**2)
        self.assertAlmostEqual(total,flux(w,10)+decay,places=12)
    def test_translation_invariance(self):
        c=coefficients(WITNESS_MODES,10);p=np.arange(6)*.65
        self.assertAlmostEqual(evaluate(c,p),evaluate(c,p+np.array(WITNESS_MODES)@np.array([.3,-.8])),places=12)
    def test_smith_torus_criterion(self):
        self.assertTrue(triangle_attainable([[2]],[-1]))
        self.assertTrue(triangle_attainable([[1,1],[1,-1]],[1,-1]))
        self.assertFalse(triangle_attainable([[1],[1]],[1,-1]))
        self.assertTrue(triangle_attainable([[1],[1]],[-1,-1]))
        self.assertFalse(triangle_attainable(WITNESS_ROWS,[1,1,-1,1]))
        self.assertTrue(triangle_attainable(WITNESS_ROWS,[1,1,1,1]))
    def test_hybrid_general_relation(self):
        modes=[(x,y) for x in range(4) for y in range(-3,4) if x>0 or (x==0 and y>0)]
        c=coefficients(modes,9);cert=pack_hybrid_certificate(c)
        self.assertEqual(cert['candidate_cycles'],0)
        self.assertGreater(cert['general_deficit'],0)
        self.assertTrue(verify_hybrid(c,cert))
        self.assertFalse(triangle_attainable(list(c),list(c.values())))
    def test_standalone_verifier(self):
        import json,importlib.util
        path=Path(__file__).resolve().parents[1]/'verify_certificates.py'
        spec=importlib.util.spec_from_file_location('standalone',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        items=json.loads((path.parent/'results/packing_certificates.json').read_text())
        for item in items:module.verify(item)
        item=items[-1];item['coefficients'][0]='1/10000000000000000'
        with self.assertRaises(AssertionError):module.verify(item)

    def test_invalid_input(self):
        with self.assertRaises(ValueError):coefficients([(0,0)],1)
        with self.assertRaises(ValueError):coefficients([(1,0),(1,0)],1)
        with self.assertRaises(ValueError):coefficients([(1,0)],1,amplitudes=[-1])

if __name__=='__main__':unittest.main(verbosity=2)
