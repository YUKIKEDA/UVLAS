"""80-digit continuum work with exact transverse polynomial moments."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from decimal import Decimal as D,localcontext
import itertools
import numpy as np
import developable_general_reference as original
def power(value,n):return D(1) if n==0 else value**n
def moment(power):return D(0) if power%2 else D(1)/(D(2)**power*(power+1))
def stress_body(kind,k,x,z,lam,mu):
    base=1+k*k*x*x;j=base-k*z;v=D(1) if kind=='quadratic' else x*x;vp=D(0) if kind=='quadratic' else 2*x;f=-z*v/j;fx=-z*(vp*j-v*2*k*k*x)/j**2;fz=-v*base/j**2
    T=[lam*base+2*mu,lam*base,lam*base+2*mu*k*k*x*x,2*mu*k*x];Tx=[lam*2*k*k*x,lam*2*k*k*x,lam*2*k*k*x+4*mu*k*k*x,2*mu*k];S=[f*t for t in T];Sx=[fx*t+f*dt for t,dt in zip(T,Tx)];Sz=[fz*t for t in T];b0=-(Sx[0]+k*x*Sx[3]-k*x*Sz[0]+(1-k*z)*Sz[3])/j;b2=-(Sx[3]+k*x*Sx[2]-k*x*Sz[3]+(1-k*z)*Sz[2])/j
    return S,[b0,D(0),b2],j,base,v
def weak(kind,k0,h0,orders):
    with localcontext() as ctx:
        ctx.prec=80;k=D(str(k0));h=D(str(h0));lam,mu=original.constants(D);p,w=original.decimal_gauss(orders[0]);xp=[(v+1)/2 for v in p];xw=[v/2 for v in w];p,w=original.decimal_gauss(orders[2]);zp=[v*h/2 for v in p];zw=[v*h/2 for v in w];indices=list(itertools.product(range(4),range(4),range(3)));internal=[[D(0)]*3 for _ in indices];volume=[[D(0)]*3 for _ in indices];outer=[[D(0)]*3 for _ in indices];reaction=[[D(0)]*3 for _ in indices];norms=[[D(0)]*3 for _ in indices];reference=D(0)
        for x,wx in zip(xp,xw):
            for z,wz in zip(zp,zw):
                weight=wx*wz;S,body,j,base,v=stress_body(kind,k,x,z,lam,mu);reference+=weight*(lam+2*mu)*z*z*v*v*base*base/j
                for index,(a,b,c) in enumerate(indices):
                    factor=D('.5') if c==2 else D(1);phi=factor*power(x,a)*power(z,c);px=factor*a*power(x,a-1)*power(z,c) if a else D(0);pz=factor*c*power(x,a)*power(z,c-1) if c else D(0);gx=(px-k*x*pz)/j;gz=(k*x*px+(1-k*z)*pz)/j;Y=moment(b);Yd=(D('.5')**b-(-D('.5'))**b);internal[index][0]+=weight*j*(gx*S[0]+gz*S[3])*Y;internal[index][2]+=weight*j*(gx*S[3]+gz*S[2])*Y;iy=weight*j*phi*S[1]*Yd;internal[index][1]+=iy;outer[index][1]+=iy
                    for alpha in [0,2]:volume[index][alpha]+=weight*j*phi*body[alpha]*Y
                    gx2=gx*gx*moment(2*b);gz2=gz*gz*moment(2*b);gy2=phi*phi*b*b*moment(2*b-2) if b else D(0);g2=[gx2,gy2,gz2]
                    for alpha in range(3):norms[index][alpha]+=weight*j*(mu*sum(g2)+(lam+mu)*g2[alpha])
        for axis in [0,2]:
            for sign in [-1,1]:
                samples=zip(zp,zw) if axis==0 else zip(xp,xw)
                for variable,weight in samples:
                    x=(D(0) if sign<0 else D(1)) if axis==0 else variable;z=variable if axis==0 else sign*h/2;S,_,_,_,_=stress_body(kind,k,x,z,lam,mu)
                    if axis==0:t=[sign*(S[0]+k*x*S[3]),D(0),sign*(S[3]+k*x*S[2])]
                    else:t=[sign*(-k*x*S[0]+(1-k*z)*S[3]),D(0),sign*(-k*x*S[3]+(1-k*z)*S[2])]
                    dest=reaction if axis==0 and sign<0 else outer
                    for index,(a,b,c) in enumerate(indices):
                        factor=D('.5') if c==2 else D(1);phi=factor*power(x,a)*power(z,c)*moment(b)
                        for alpha in [0,2]:dest[index][alpha]+=weight*phi*t[alpha]
        error=[];scale=D(0);maxc=D(0);allowed=D(0)
        for index,(a,b,c) in enumerate(indices):
            for alpha in range(3):
                r=internal[index][alpha]-volume[index][alpha]-outer[index][alpha]-reaction[index][alpha];error.append(r*r);scale+=internal[index][alpha]**2
                if norms[index][alpha]>0:maxc=max(maxc,abs(r)/(reference*norms[index][alpha]).sqrt())
                if a>0:allowed=max(allowed,abs(reaction[index][alpha]))
        return dict(internal=np.array([[float(v) for v in row] for row in internal]),relative=float((sum(error)/scale).sqrt()),cauchy=float(maxc),allowed_reaction=float(allowed),reference=float(reference))
