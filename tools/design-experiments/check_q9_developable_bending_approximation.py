"""Full-energy best approximation on an exact Q9 developable reference geometry."""
from pathlib import Path
import itertools,json,time
import numpy as np
import q9_developable_bending_approximation as method
import q9_developable_bending_bounds as certificate
core=method.core;ROOT=core.ROOT;F=core.F
def gate(cases,route):
    selected=[c for c in cases if c['kappa']==.4 and c['amplitude']==.2 and c['n']==8];finest=[]
    for c in selected:
        x=next(x for x in c['solutions'] if x['route']==route and x['factor']==2);finest.append(dict(h=c['h'],lower=x['bounds']['original9_lower_relative'],upper=x['observation11']['error_upper']))
    coarse=[c for c in cases if c['kappa']==.4 and c['amplitude']==.2 and c['n']==4];a=next((c for c in coarse if c['h']==.01),None);b=next((c for c in coarse if c['h']==.0001),None);ratio=None
    if a and b:
        x=next(x for x in a['solutions'] if x['route']==route and x['factor']==2);y=next(x for x in b['solutions'] if x['route']==route and x['factor']==2);lo0=x['bounds']['original9_lower_relative'];up0=x['observation11']['error_upper'];lo1=y['bounds']['original9_lower_relative'];up1=y['observation11']['error_upper'];ratio=dict(lower=lo1/up0 if up0>1e-12 else None,upper=up1/lo0 if lo0>1e-12 else None)
    failed=any(x['lower']>.01 for x in finest) or bool(ratio and ratio['lower'] is not None and ratio['lower']>2)
    passed=len(finest)==2 and all(x['upper']<=.01 for x in finest) and ratio is not None and ratio['upper'] is not None and ratio['upper']<=2
    return dict(status='failed' if failed else ('passed' if passed else 'unresolved'),finest=finest,thin_amplification=ratio)
def main():
    target=F/'q9-developable-bending-approximation-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-developable-bending-reference-result.json';pr=json.loads(parent.read_text());assert pr['all_criteria_pass'];ref={(c['kappa'],c['h']):float(c['decimal_energy']) for c in pr['cases']};test=certificate.synthetic();test.pop('feasible_r');test.pop('feasible_p');start=time.monotonic();cases=[];stop=None;out=ROOT/'build/analysis'/('q9-developable-bending-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    configs=list(itertools.product([0.,.4],[0.,.2],[1,2,4],[.01,.001,.0001]))+[(.4,.2,8,h) for h in [.001,.0001]]
    try:
        for k,amplitude,n,h in configs:
            states={};a=method.assemble(k,n,h,amplitude,states);Rsave={};records=[];fields=[];saved_r=[];saved_p=[];globals_={};spectra={}
            for route in ['raw','cg']:
                factors=[2.,4.] if k==.4 and n==4 and h==.0001 else [2.];R,qrlog,solutions,spectrum=method.solve_route(a,route,factors);Rsave[route]=R;globals_[route]=qrlog;spectra[route]=dict(position_condition=spectrum['position_condition'],singular_min=float(spectrum['direction_singular_values'][-1]),singular_max=float(spectrum['direction_singular_values'][0]));local=[x for x in a['locals'] if x['route']==route];delta=float(np.sqrt(sum(x['qr']['error_upper']**2 for x in local))/np.sqrt(a['reference_norm']));eta=max(x['qr']['orthogonality_upper'] for x in local);zero_upper=1/np.sqrt(1-method.rounding.gamma(6*729*len(a['cells'])+10))
                for fit in solutions:
                    r,p=fit.pop('r'),fit.pop('p');bounds=certificate.certify(R,a['routes'][route]['W'].shape[1],a['M'],r,p,fit['h1']*0+fit['factor']*fit['reference_h1'],qrlog,delta,eta,zero_upper);fr=bounds.pop('feasible_r');fp=bounds.pop('feasible_p');u=method.full_field(a,fr,fp);fields.append((route,u));saved_r.append(fr);saved_p.append(fp);fit.update(route=route,bounds=bounds,position_error9=method.position_error(k,a,u,9),position_error11=method.position_error(k,a,u,11));records.append(fit)
            # Full rows are no longer needed after both global QR solves.
            a['routes']={};reference9,obs9=method.observe(k,h,a,states,fields,9);reference11,obs11=method.observe(k,h,a,states,fields,11)
            for fit,o9,o11 in zip(records,obs9,obs11):fit['observation9']=o9;fit['observation11']=o11;fit['compression_difference']=abs(o9['error']-fit['bounds']['qr_upper_relative']);fit['integration_difference']=abs(o11['error']-o9['error'])
            reference_error=max(abs(reference9/(2*ref[k,h])-1),abs(reference11/(2*ref[k,h])-1));checks=dict(geometry=max(a['geometry_position'],a['geometry_jacobian'])<=1e-11,edge=a['edge_valid'],local_rank=all(x['ranks']==[75,75,75] for x in a['locals']),rigid=max(x['rigid'] for x in a['locals'])<=1e-11,kkt=all(max(x['kkt'],x['kkt_r'],x['kkt_p'])<=1e-9 for x in records),sphere=all(x['constraint_error']<=1e-10 and x['bounds']['feasible'] for x in records),bounds=all(x['bounds']['objective_order'] for x in records),reference=reference_error<=1e-10,integration=all(x['integration_difference']<=1e-7 for x in records),compression=all(x['compression_difference']<=1e-7 for x in records),position_observation=all(abs(x['position_error9']-x['position_error11'])<=1e-7 for x in records),flat=k!=0 or all(x['observation11']['error_upper']<=1e-8 for x in records))
            path=out/f'k{k:g}-a{amplitude:g}-n{n}-h{h:g}.npz';np.savez_compressed(path,Rraw=Rsave['raw'],Rcg=Rsave['cg'],M=a['M'],xy=a['xy'],coords=a['coords'],reference=a['reference'],free=a['free'],r=np.array(saved_r),p=np.array(saved_p));record=dict(kappa=k,amplitude=amplitude,n=n,h=h,geometry_position=a['geometry_position'],geometry_jacobian=a['geometry_jacobian'],edge_condition=a['edge_condition'],reference_error=reference_error,reference_energy=reference11/2,reference_interpolation=a['interpolation'],interpolation_position_error=method.position_error(k,a,a['reference'],11),reduction_change=a['reduction_change'],local_diagnostics=a['locals'],global_qr=globals_,spectra=spectra,solutions=records,checks={key:bool(value) for key,value in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=core.digest(path));cases.append(record)
            print(json.dumps(dict(kappa=k,amplitude=amplitude,n=n,h=h,errors={x['route']+str(x['factor']):[x['bounds']['original9_lower_relative'],x['observation11']['error_upper']] for x in records},checks=record['checks'],seconds=time.monotonic()-start)),flush=True)
            del states,a,Rsave,fields
            if time.monotonic()-start>1200 or core.memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
        assert sum(p.stat().st_size for p in out.iterdir())<=256*1024**2
    except Exception as exc:stop=repr(exc)
    complete=stop is None and len(cases)==38 and sum(len(c['solutions']) for c in cases)==80
    result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),protocol_sha256=core.digest(F/'q9-developable-bending-approximation-experiment.md'),definition_sha256=core.digest(ROOT/'docs/q9-developable-bending-reference.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [method,certificate,core,core.method,core.method.base,method.material,core.geometry]},cases=cases,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),numeric_checks_pass=bool(complete and all(all(c['checks'].values()) for c in cases)),candidate_gate=gate(cases,'cg'),raw_gate=gate(cases,'raw'),product_shell_adopted=False,load_equilibrium_reference_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({key:value for key,value in result.items() if key not in ['cases','dependencies']}))
    if not result['numeric_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
