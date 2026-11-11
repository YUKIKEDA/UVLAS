"""MR-R1: exact rational modal integration plus independent Decimal LDL solve."""
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal as D,localcontext
import hashlib,json,math,time
import check_q9_independent_load as material
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fq(x):return Q.from_float(float(x))
def dec(x):return D(x.numerator)/D(x.denominator)
def integrate(h,p,C):
    G=[[{} for _ in range(4)] for _ in range(6)]
    G[0][0]={1:-p};G[1][1]={1:-p};G[2][3]={1:Q(1)}
    G[3][1]={0:Q(1)};G[3][2]={0:p};G[3][3]={2:p/2}
    G[4][0]={0:Q(1)};G[4][2]={0:p};G[4][3]={2:p/2}
    G[5][0]={1:p};G[5][1]={1:p}
    K=[[Q(0) for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            for a in range(6):
                for b in range(6):
                    if C[a][b]:
                        for degree,x in G[a][i].items():
                            for power,y in G[b][j].items():
                                k=degree+power
                                if k%2==0:K[i][j]+=x*y*C[a][b]*h**(k+1)/Q(4*2**k*(k+1))
    return K
def closed(h,p,d,l,u):
    i0=h;i2=h**3/12;i4=h**5/80
    a=u*i0+(d+u)*p*p*i2;b=(l+u)*p*p*i2;c=u*p*i0;e=(u/2-l)*p*i2
    return [[v/4 for v in row] for row in [[a,b,c,e],[b,a,c,e],[c,c,2*u*p*p*i0,u*p*p*i2],[e,e,u*p*p*i2,d*i2+u*p*p*i4/2]]]
def rational_solve(K,b):
    a=[list(row)+[v] for row,v in zip(K,b)];n=len(b)
    for i in range(n):
        k=max(range(i,n),key=lambda j:abs(a[j][i]));a[i],a[k]=a[k],a[i];pivot=a[i][i];assert pivot
        a[i]=[v/pivot for v in a[i]]
        for j in range(n):
            if j!=i:
                factor=a[j][i];a[j]=[v-factor*w for v,w in zip(a[j],a[i])]
    return [a[i][-1] for i in range(n)]
def ldl(K,b):
    n=len(b);L=[[D(0) for _ in range(n)] for _ in range(n)];diag=[]
    for i in range(n):
        L[i][i]=D(1);v=K[i][i]-sum((L[i][j]**2*diag[j] for j in range(i)),D(0));assert v>0;diag.append(v)
        for k in range(i+1,n):L[k][i]=(K[k][i]-sum((L[k][j]*L[i][j]*diag[j] for j in range(i)),D(0)))/v
    y=[]
    for i in range(n):y.append(b[i]-sum((L[i][j]*y[j] for j in range(i)),D(0)))
    z=[v/d for v,d in zip(y,diag)];x=[D(0)]*n
    for i in reversed(range(n)):x[i]=z[i]-sum((L[j][i]*x[j] for j in range(i+1,n)),D(0))
    return x
def pi_rational():
    def atan(x):
        result=Q(0);k=0
        while True:
            term=x**(2*k+1)/Q(2*k+1);result+=(-1)**k*term;k+=1
            tail=x**(2*k+1)/Q(2*k+1)
            if tail<Q(1,10**105):return result,tail
    a,ta=atan(Q(1,5));b,tb=atan(Q(1,239))
    return 16*a-4*b,16*ta+4*tb
def main():
    output=F/'plate-modal-exact-reference-result.json'
    if output.exists():raise FileExistsError(output)
    start=time.monotonic();parent=F/'q9-distorted-plate-result.json';data=json.loads(parent.read_text());C=[[fq(v) for v in row] for row in material.elastic];pi=fq(math.pi);records=[];piq,pitail=pi_rational()
    for h64 in [.01,.001]:
        old=[c['reference'] for c in data['cases'] if c['h']==h64];assert len(old)==4 and all(r==old[0] for r in old)
        h=fq(h64);pressure=fq(old[0]['pressure']);K=integrate(h,pi,C);direct=closed(h,pi,C[0][0],C[0][1],C[3][3]);assert K==direct
        b=[Q(0),Q(0),pressure/4,pressure*h*h/32];x=rational_solve(K,b);J=sum((v*w for v,w in zip(b,x)),Q(0))
        residual=[sum((a*v for a,v in zip(row,x)),Q(0))-rhs for row,rhs in zip(K,b)];assert all(v==0 for v in residual)
        decimal_paths={};reference=[float(v) for v in x];rounded_J=float(J)
        for digits,limit in [(80,D('1e-50')),(100,D('1e-65'))]:
            with localcontext() as ctx:
                ctx.prec=digits;kd=closed(dec(h),dec(pi),dec(C[0][0]),dec(C[0][1]),dec(C[3][3]));bd=[D(0),D(0),dec(pressure)/4,dec(pressure)*dec(h)**2/32];xd=ldl(kd,bd);jd=sum((v*w for v,w in zip(bd,xd)),D(0))
            with localcontext() as ctx:
                ctx.prec=150;errors=[abs(v/dec(w)-1) for v,w in zip(xd,x)]+[abs(jd/dec(J)-1)];assert max(errors)<=limit
            assert [float(v).hex() for v in xd]==[v.hex() for v in reference] and float(jd).hex()==rounded_J.hex()
            decimal_paths[str(digits)]=dict(coefficients=list(map(str,xd)),compliance=str(jd),maximum_relative_difference=str(max(errors)))
        rounding=[float(abs(v-fq(float(v)))/fq(math.ulp(float(v)))) for v in x+[J]];assert max(rounding)<=.5
        with localcontext() as ctx:
            ctx.prec=100;kh=closed(dec(h),dec(piq),dec(C[0][0]),dec(C[0][1]),dec(C[3][3]));bd=[D(0),D(0),dec(pressure)/4,dec(pressure)*dec(h)**2/32];xp=ldl(kh,bd)
            pi_sensitivity=float(abs(xp[2]/dec(x[2])-1));pi_truncation_bound=str(dec(pitail))
        records.append(dict(h=h64,inputs_hex=dict(h=h64.hex(),pi=math.pi.hex(),pressure=old[0]['pressure'].hex(),elastic=[[float(v).hex() for v in row] for row in C]),
            matrix_rational=[[str(v) for v in row] for row in K],rhs_rational=list(map(str,b)),coefficients_rational=list(map(str,x)),compliance_rational=str(J),
            rounded_reference=dict(coefficients=reference,center=reference[2],pressure=old[0]['pressure'],compliance=rounded_J),decimal_paths=decimal_paths,rounding_errors_in_ULP=rounding,
            exact_residual_zero=True,closed_matrix_matches=True,pi_series_truncation_bound=pi_truncation_bound,pi_input_center_sensitivity=pi_sensitivity,
            old_center_relative_difference=abs(old[0]['center']/reference[2]-1)))
        assert time.monotonic()-start<=30 and memory.memory_bytes()<=512*1024**2
    result=dict(version='MR-R1',source_sha256=digest(__file__),protocol_sha256=digest(ROOT/'docs/plate-modal-reference-accuracy-contract.md'),parent_sha256=digest(parent),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [material,memory]},cases=records,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=True,product_shell_adopted=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=True,seconds=result['seconds'],cases=[dict(h=c['h'],old_center_difference=c['old_center_relative_difference'],pi_input_sensitivity=c['pi_input_center_sensitivity'],decimal100_difference=c['decimal_paths']['100']['maximum_relative_difference']) for c in records]),indent=2))
if __name__=='__main__':main()
