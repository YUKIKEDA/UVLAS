"""Design experiment only: immutable step fields and Decimal geometric overlay."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from dataclasses import dataclass
from decimal import Decimal,localcontext
from bisect import bisect_right
import math
import numpy as np

def D(x):return x if isinstance(x,Decimal) else Decimal.from_float(float(x))
def line(a,b):
    with localcontext() as ctx:
        ctx.prec=80;s=(D(b[0])-D(a[0]))/(D(b[1])-D(a[1]))
        return (D(a[0])-s*D(a[1]),s)
def at(l,y):return l[0]+l[1]*y

@dataclass(frozen=True)
class Domain:
    key:str
    pieces:tuple # lo, hi, LE affine pair, TE affine pair
    def piece(self,y):return next(s for s in self.pieces if s[0]<=y<=s[1])
    @property
    def low(self):return self.pieces[0][0]
    @property
    def high(self):return self.pieces[-1][1]

@dataclass(frozen=True)
class Strip:
    low:Decimal
    high:Decimal
    fronts:tuple
    strengths:tuple

@dataclass(frozen=True)
class Field:
    domain:Domain
    strips:tuple
    basis:str='horseshoe'
    def __post_init__(self):
        with localcontext() as ctx:
            ctx.prec=80
            if self.basis!='horseshoe':raise ValueError('basis_mismatch')
            if not self.strips or self.strips[0].low!=self.domain.low or self.strips[-1].high!=self.domain.high:raise ValueError('incomplete_domain')
            previous=self.domain.low
            for s in self.strips:
                if s.low!=previous or s.high<=s.low:raise ValueError('gap_overlap_or_degenerate_strip')
                previous=s.high
                if len(s.fronts)!=len(s.strengths) or not s.fronts:raise ValueError('shape')
                if not all(math.isfinite(float(v)) for v in s.strengths):raise ValueError('nonfinite')
                if not all(v.is_finite() for l in s.fronts for v in l):raise ValueError('nonfinite')
                for y in [s.low,s.high]:
                    xs=[at(l,y) for l in s.fronts]
                    if any(b<=a for a,b in zip(xs[:-1],xs[1:])):raise ValueError('front_order_or_crossing')
                    piece=self.domain.piece(y);le=at(piece[2],y);te=at(piece[3],y)
                    if xs[0]<le-D('1e-11') or xs[-1]>te+D('1e-11'):raise ValueError('outside_wing')
    def active(self,y):return self.strips[bisect_right([s.high for s in self.strips],y)]

def wing_domain(case):
    with localcontext() as ctx:
        ctx.prec=80
        t=D(math.tan(math.radians(35)));s=D(.225);c=D(1.5)
        return Domain(case,((D(-4),D(0),(D(0),-t),(c,-t+s)),(D(0),D(4),(D(0),t),(c,t-s))))

def rectangle():return Domain('manufactured',((D(0),D(1),(D(0),D(0)),(D(1),D(0))),))

def integral_exact(field):
    terms=[]
    with localcontext() as ctx:
        ctx.prec=80
        for s in field.strips:
            ys=sorted(set([s.low,s.high]+[p[0] for p in field.domain.pieces if s.low<p[0]<s.high]))
            for lo,hi in zip(ys[:-1],ys[1:]):
                te=field.domain.piece((lo+hi)/2)[3]
                for l,g in zip(s.fronts,s.strengths):
                    amount=(at(te,lo)-at(l,lo)+at(te,hi)-at(l,hi))*(hi-lo)/2
                    terms.append(float(amount)*g)
    return math.fsum(terms)

def overlay(fields,observe=False):
    if not fields or any(f.domain!=fields[0].domain for f in fields):raise ValueError('domain_mismatch')
    domain=fields[0].domain;areas=[];values=[];cells=[];crossings=0;rounded=Decimal(0)
    with localcontext() as ctx:
        ctx.prec=80
        ybounds=[domain.low+(domain.high-domain.low)*D(j)/D(20) for j in range(21)]
        ys={domain.low,domain.high}
        for p in domain.pieces:ys.update(p[:2])
        for field in fields:
            for s in field.strips:ys.update([s.low,s.high])
        if observe:ys.update(ybounds)
        ys=sorted(ys)
        for lo,hi in zip(ys[:-1],ys[1:]):
            mid=(lo+hi)/2;piece=domain.piece(mid);le,te=piece[2:]
            lines=[le,te];labels=[(-1,0.),(-2,0.)]
            for k,f in enumerate(fields):
                s=f.active(mid);lines.extend(s.fronts);labels.extend((k,g) for g in s.strengths)
            if observe:
                for j in range(1,16):
                    q=D(j)/D(16);lines.append(tuple(a+q*(b-a) for a,b in zip(le,te)));labels.append((-3,j))
            a=np.array([float(l[0]) for l in lines]);b=np.array([float(l[1]) for l in lines])
            ii,jj=np.triu_indices(len(lines),1);da=a[ii]-a[jj];db=b[ii]-b[jj]
            roots=np.divide(-da,db,out=np.full(len(da),np.inf),where=db!=0)
            near=(abs(da)<=1e-12)&(abs(db)<=1e-12)
            possible=((roots>=float(lo)-1e-12)&(roots<=float(hi)+1e-12))|near
            cuts={lo,hi}
            for i,j in zip(ii[possible],jj[possible]):
                denominator=lines[i][1]-lines[j][1]
                if denominator:
                    root=-(lines[i][0]-lines[j][0])/denominator
                    if lo<root<hi:cuts.add(root)
            cuts=sorted(cuts);crossings+=len(cuts)-2
            for y0,y1 in zip(cuts[:-1],cuts[1:]):
                ym=(y0+y1)/2;order=sorted(range(len(lines)),key=lambda j:at(lines[j],ym))
                contributions=[[] for _ in fields];current=[0.]*len(fields);inside=False;chord=0
                span=min(bisect_right(ybounds,ym)-1,19)
                for left,right in zip(order[:-1],order[1:]):
                    label,g=labels[left]
                    if label==-1:inside=True
                    elif label==-2:inside=False
                    elif label==-3:chord=int(g)
                    else:
                        contributions[label].append(g);current[label]=math.fsum(contributions[label])
                    if not inside:continue
                    width0=at(lines[right],y0)-at(lines[left],y0)
                    width1=at(lines[right],y1)-at(lines[left],y1)
                    area=(width0+width1)*(y1-y0)/2
                    if area<0:
                        if abs(area)>Decimal('1e-60'):raise ArithmeticError('negative_area')
                        rounded+=abs(area);continue
                    if area==0:continue
                    areas.append(float(area));values.append(tuple(current))
                    if observe:cells.append(span*16+chord)
    a=np.array(areas);v=np.array(values)
    if not len(a) or not np.all(a>0):raise ArithmeticError('empty_or_underflow_area')
    integrals=[math.fsum(a*v[:,j]) for j in range(len(fields))]
    result=dict(area=math.fsum(a),regions=len(a),crossings=crossings,rounded_negative_area=str(rounded),
        minimum_region_area=float(np.min(a)),regions_below_1e60=int(np.count_nonzero(a<1e-60)),
        integrals=integrals,integral_preservation=max(abs(x-integral_exact(f)) for x,f in zip(integrals,fields)),
        norms=[metrics(a,v[:,j]) for j in range(len(fields))])
    if observe:
        cell=np.array(cells);cell_area=[math.fsum(a[cell==j]) for j in range(320)]
        if min(cell_area)<=0:raise ArithmeticError('uncovered_observation_cell')
        result['cell_areas']=cell_area
        result['cell_integrals']=[[math.fsum(a[cell==j]*v[cell==j,k]) for j in range(320)] for k in range(len(fields))]
    return result,a,v

def metrics(a,v):return dict(l1=math.fsum(a*abs(v)),l2=math.sqrt(math.fsum(a*v*v)),linf=float(np.max(abs(v))))

def diagnostic(differences,cosines,scale):
    d=np.array(differences,dtype=float)
    out=dict(differences=d.tolist(),difference_cosines=cosines,status='not_asymptotic',orders=None)
    if max(d)<=1e-10*scale:
        out['status']='roundoff_candidate';return out
    if min(d)<=0:return out
    p=np.log2(d[:-1]/d[1:]);out['orders']=p.tolist()
    if all(x is not None and x>=.9 for x in cosines) and np.all(np.diff(d)<0) and .5<=p[-1]<=6 and abs(p[-1]-p[-2])<=.5:out['status']='asymptotic_candidate'
    return out
