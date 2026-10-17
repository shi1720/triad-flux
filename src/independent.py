"""Independent physical-space and direct-convolution checks (no flux grouping)."""
import numpy as np

def state(modes, amplitudes, phases, n):
    w=np.zeros((n,n),complex)
    for k,a,p in zip(modes,amplitudes,phases):
        w[k[0]%n,k[1]%n]=a*np.exp(1j*p)
        w[-k[0]%n,-k[1]%n]=a*np.exp(-1j*p)
    return w

def rhs(w, viscosity=0., dealias=True):
    n=w.shape[0];k=np.rint(np.fft.fftfreq(n)*n);kx,ky=np.meshgrid(k,k,indexing='ij')
    r=kx*kx+ky*ky;safe=np.where(r==0,1,r)
    phys=lambda a: np.fft.ifft2(a*n*n).real
    ux=phys(1j*ky*w/safe);uy=phys(-1j*kx*w/safe)
    wx=phys(1j*kx*w);wy=phys(1j*ky*w)
    nonlinear=-np.fft.fft2(ux*wx+uy*wy)/(n*n)
    if dealias: nonlinear*=((abs(kx)<n/3)&(abs(ky)<n/3))
    return nonlinear-viscosity*r*w

def flux(w,cutoff2,observable='energy'):
    n=w.shape[0];k=np.rint(np.fft.fftfreq(n)*n);kx,ky=np.meshgrid(k,k,indexing='ij')
    r=kx*kx+ky*ky;mask=(r>0)&(r<=cutoff2)
    weight=1/np.where(r==0,1,r) if observable=='energy' else np.ones_like(r)
    return -float(np.sum(weight[mask]*(w.conj()*rhs(w))[mask].real))

def direct_flux(modes,amps,phases,cutoff2,observable='energy'):
    vals={}
    for p,a,phi in zip(modes,amps,phases):
        vals[tuple(p)]=a*np.exp(1j*phi);vals[tuple(-x for x in p)]=a*np.exp(-1j*phi)
    total=0.
    for k,wk in vals.items():
        kn=sum(x*x for x in k)
        if kn>cutoff2:continue
        dw=0j
        for p,wp in vals.items():
            q=(k[0]-p[0],k[1]-p[1])
            if q not in vals:continue
            dw-=(p[0]*q[1]-p[1]*q[0])/sum(x*x for x in p)*wp*vals[q]
        total-=(wk.conjugate()*dw).real*(1/kn if observable=='energy' else 1)
    return float(total)

def rk4(w,dt,viscosity):
    a=rhs(w,viscosity);b=rhs(w+.5*dt*a,viscosity)
    c=rhs(w+.5*dt*b,viscosity);d=rhs(w+dt*c,viscosity)
    return w+dt*(a+2*b+2*c+d)/6
