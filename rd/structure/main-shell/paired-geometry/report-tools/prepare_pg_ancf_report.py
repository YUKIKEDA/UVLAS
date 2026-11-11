"""Read frozen evidence; generate report tables, figures and independent identities."""
from pathlib import Path
import argparse,csv,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'build/paper-tools'))
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=ROOT/'docs/papers/pg-ancf-assets'
NAMES=['mass-paired-quality','mass-paired-coefficient-input','mass-paired-reference',
       'paired-lift-exact','mass-paired-body','mass-paired-interface','mass-paired-mechanics',
       'mass-paired-orientation','paired-shell-r5-domain-v2']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',help='New root-relative directory for independent regeneration')
    args=parser.parse_args()
    if args.output_dir:
        OUT=(ROOT/args.output_dir).resolve()
        if not OUT.is_relative_to(ROOT.resolve()) or OUT==ROOT.resolve():
            raise SystemExit('Output must stay below the workspace root')
        if OUT.exists() and any(OUT.iterdir()):
            raise SystemExit('Regeneration requires a new or empty output directory')
    elif (ROOT/'docs/papers/source-hashes.json').exists():
        raise SystemExit('Report v1 is frozen; use --output-dir with a new directory')
    start=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True)
    data={};sources={}
    for name in NAMES:
        p=ROOT/f'rd/{name}-result.json';d=json.loads(p.read_text());assert d['passed'];data[name]=d;sources[p.relative_to(ROOT).as_posix()]=sha(p)
        for rel,h in d.get('sources',{}).items():assert sha(ROOT/rel)==h,rel;sources[rel]=h
        for row in d.get('records',[]):
            if 'artifact' in row:assert sha(ROOT/row['artifact'])==row['artifact_sha256']
    s,t,u=sp.symbols('s t u',real=True);a=sp.symbols('a',real=True);w=sp.symbols('w',positive=True)
    f=sp.Function('f');k=sp.Function('k');I0=sp.Integral(k(u),(u,0,s));I1=sp.Integral(u*k(u),(u,0,s));b=s*I0-I1
    J=sp.Integral(f(u,t)*k(u),(u,0,s));Jy=sp.Integral(sp.diff(f(u,t),t)*k(u),(u,0,s));Jsy=sp.Integral(u*sp.diff(f(u,t),t)*k(u),(u,0,s))
    v=sp.Matrix([(-f(s,t)*I0+J)/w,-s*Jy+Jsy-a*(-f(s,t)*I0+J)/w,b]);r=sp.Matrix([w*s+a*t,t,f(s,t)]);As=r.diff(s);At=r.diff(t)
    identities={}
    for name,value in [('membrane_ss',As.dot(v.diff(s))),('membrane_st',(As.dot(v.diff(t))+At.dot(v.diff(s)))/2),('membrane_tt_retained',At.dot(v.diff(t))-sp.diff(-s*Jy+Jsy,t))]:
        residual=sp.simplify(value);assert residual==0,(name,residual);identities[name]=str(residual)
    p,q,L1,L2=sp.symbols('p q L1 L2',real=True);g=1+p*p+q*q
    normal=sp.Matrix([-p,-q,1])/sp.sqrt(g);director=-sp.Matrix([(1+q*q)*L1-p*q*L2,(1+p*p)*L2-p*q*L1,p*L1+q*L2])/g**sp.Rational(3,2)
    for name,value in [('director_tangent_1',sp.Matrix([1,0,p]).dot(director)+L1/sp.sqrt(g)),('director_tangent_2',sp.Matrix([0,1,q]).dot(director)+L2/sp.sqrt(g)),('director_normal',normal.dot(director))]:
        residual=sp.simplify(value);assert residual==0;identities[name]=str(residual)
    # A zero director cannot cancel a nonzero transverse slope even on a flat graph.
    neg=sp.diff(s*s/2,s).subs(s,sp.Rational(1,2))/2;assert neg!=0
    symbolic=dict(passed=True,identities=identities,zero_director_negative_control=str(neg),scope='independent symbolic differential/algebraic identities; not uniform convergence or novelty')
    (OUT/'symbolic-checks.json').write_text(json.dumps(symbolic,indent=2),encoding='utf-8')
    rows=[];geometries=[]
    for record in data['mass-paired-quality']['records']:
        key=tuple(record['parameters'])
        if key not in geometries:geometries.append(key)
        obs=record['observation'];rows.append(dict(group='G'+str(geometries.index(key)+1),family=record['family'],n=record['n'],h=record['h'],dofs=3*record['ns'],energy=max(obs['error']),position=max(obs['position']),stationarity=max(obs['stationarity']),seconds=record['cold_case_seconds']))
    for record in data['mass-paired-coefficient-input']['records']:
        obs=record['observation'];rows.append(dict(group='G7',family='coefficient-map',n=record['n'],h=record['h'],dofs=record['stiffness_ranks'][0],energy=max(obs['error']),position=max(obs['position']),stationarity=max(obs['stationarity']),seconds=record['seconds']))
    assert len(rows)==28 and len(geometries)==6
    with (OUT/'quality-cases.csv').open('w',encoding='utf-8',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    amplification=[]
    for gid in ['G'+str(j) for j in range(1,8)]:
        original=[r for r in rows if r['group']==gid];finest=max(r['n'] for r in original);fine=[r for r in original if r['n']==finest]
        # Preserve original per-RHS thin ratios, rather than dividing two maxima.
        doc=data['mass-paired-coefficient-input'] if gid=='G7' else data['mass-paired-quality']
        comps=[x for x in doc['comparisons'] if x['kind']=='thin_amplification' and (gid=='G7' or tuple(x['parameters'])==geometries[int(gid[1:])-1])]
        ratio=max(max(x['values']) for x in comps)
        amplification.append(dict(group=gid,n=finest,max_energy=max(r['energy'] for r in fine),max_position=max(r['position'] for r in fine),max_thin_amplification=ratio,max_dofs=max(r['dofs'] for r in fine)))
    physical=[]
    for key,name in [('rectangle','mass-paired-body'),('taper','paired-shell-r5-domain-v2')]:
        for r0 in data[name]['records']:
            obs=r0['observation'] if key=='rectangle' else r0['physical'];physical.append(dict(domain=key,n=r0['n'],h=r0['h'],energy=max(obs['error']),position=max(obs['position']),stationarity=max(obs['stationarity']),discrete_load_agreement=r0['physical_vs_projection'] if key=='rectangle' else r0['physical_vs_weak']))
    fine=[r for r in rows if r['n']==max(q['n'] for q in rows if q['group']==r['group'])]
    summary=dict(quality_count=28,quality_rhs_count=84,quality_finest_max_energy=max(r['energy'] for r in fine),quality_finest_max_position=max(r['position'] for r in fine),quality_max_stationarity=max(r['stationarity'] for r in rows),quality_max_case_seconds=max(r['seconds'] for r in rows),groups=amplification,physical=physical,geometry_parameters=[list(x) for x in geometries],sources=sources)
    for name,d in [('summary.json',summary)]: (OUT/name).write_text(json.dumps(d,indent=2),encoding='utf-8')
    tables=['| Group | finest n | DOFs | max energy | max position | max thin ratio* |','|---|---:|---:|---:|---:|---:|']
    for r0 in amplification:tables.append(f"| {r0['group']} | {r0['n']} | {r0['max_dofs']} | {r0['max_energy']:.6e} | {r0['max_position']:.6e} | {r0['max_thin_amplification']:.6f} |")
    tables+=['','*Thin ratio is the maximum original per-RHS energy/position ratio over BOTH saved resolutions, with the original 1e-12 floor. It is not an asymptotic convergence bound.','','| Domain | n | h | energy | position | original stationarity | load agreement |','|---|---:|---:|---:|---:|---:|---:|']
    for r0 in physical:tables.append(f"| {r0['domain']} | {r0['n']} | {r0['h']:.0e} | {r0['energy']:.6e} | {r0['position']:.6e} | {r0['stationarity']:.3e} | {r0['discrete_load_agreement']:.3e} |")
    (OUT/'result-tables.md').write_text('\n'.join(tables)+'\n',encoding='utf-8')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.grid':True,'grid.alpha':.25,'savefig.bbox':'tight'})
    fig,axes=plt.subplots(1,2,figsize=(8.2,3.15));colors=plt.cm.tab10.colors
    for idx,gid in enumerate(['G'+str(j) for j in range(1,8)]):
        ax=axes[0 if idx<4 else 1]
        for h,style in [(.01,'-o'),(.0001,'--s')]:
            selected=sorted([r for r in rows if r['group']==gid and r['h']==h],key=lambda r:r['dofs']);ax.loglog([r['dofs'] for r in selected],[r['energy'] for r in selected],style,color=colors[idx],label=gid+(' (h=1e-2)' if h==.01 else ' (h=1e-4)'),markersize=4)
    for ax in axes:ax.set_xlabel('Constrained Cartesian DOFs');ax.set_ylabel('Max. relative energy error (3 RHS)');ax.legend(fontsize=6,ncol=2)
    axes[0].set_title('Quadratic reference graphs');axes[1].set_title('Nonquadratic / general coefficient graphs');fig.tight_layout();fig.savefig(OUT/'quality-refinement.svg');fig.savefig(OUT/'quality-refinement.png',dpi=220);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.8,2.7));xx=np.arange(len(amplification));ax.bar(xx,[r['max_thin_amplification'] for r in amplification],color='#286178');ax.axhline(2,color='#a64935',linestyle='--',label='Internal acceptance limit = 2');ax.set_xticks(xx,[r['group'] for r in amplification]);ax.set_ylabel('Max. thin amplification');ax.set_ylim(0,2.25);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'thin-amplification.svg');fig.savefig(OUT/'thin-amplification.png',dpi=220);plt.close(fig)
    # Independent PDE, all two thicknesses and both domains; no fitted trend line.
    fig,axes=plt.subplots(1,2,figsize=(8.2,2.7))
    for domain,color in [('rectangle','#286178'),('taper','#a64935')]:
        for h,style in [(.01,'-o'),(.0001,'--s')]:
            rr=sorted([r for r in physical if r['domain']==domain and r['h']==h],key=lambda r:r['n'])
            for ax,key in zip(axes,['energy','position']):ax.semilogy([r['n'] for r in rr],[r[key] for r in rr],style,color=color,label=f'{domain}, h={h:g}',markersize=4);ax.set_xticks([1,2]);ax.set_xlabel('Internal patch resolution n')
    axes[0].set_ylabel('Relative energy error');axes[1].set_ylabel('Relative displacement error');axes[0].legend(fontsize=6);fig.tight_layout();fig.savefig(OUT/'independent-pde.svg');fig.savefig(OUT/'independent-pde.png',dpi=220);plt.close(fig)
    manifest=dict(passed=True,meaning='report aggregation and independent identity checks only',sources=sources,script_sha256=sha(Path(__file__)),symbolic=symbolic,seconds=time.perf_counter()-start,artifacts={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='preparation.json'})
    (OUT/'preparation.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,quality_cases=len(rows),finest_max_energy=summary['quality_finest_max_energy'],symbolic_identities=len(identities),seconds=manifest['seconds'])))
if __name__=='__main__':main()
