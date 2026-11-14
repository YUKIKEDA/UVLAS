"""Global material aliases and oriented TE contract; design witness, not public API."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
def compile_bindings(patches,te_edges):
    keys=tuple(dict.fromkeys(k for p in patches for k in p['dofs']));cols={k:i for i,k in enumerate(keys)};vertices={};rows=[];ids=[];faces=[]
    for p in patches:
        N=np.asarray(p['N'],float)
        if N.shape!=(len(p['vertices']),len(p['dofs'])) or not np.isfinite(N).all() or len(set(p['dofs']))!=len(p['dofs']):raise ValueError('invalid_global_binding')
        wide=np.zeros((len(N),len(keys)))
        for j,key in enumerate(p['dofs']):wide[:,cols[key]]=N[:,j]
        for j,key in enumerate(p['vertices']):
            if key in vertices:
                if not np.array_equal(rows[vertices[key]],wide[j]):raise ValueError('alias_position_field_mismatch')
            else:vertices[key]=len(rows);rows.append(wide[j]);ids.append(key)
        for face in p['faces']:
            f=tuple(vertices[p['vertices'][j]] for j in face)
            if len(set(f))!=len(f):raise ValueError('degenerate_topological_face')
            if any(set(f)==set(old) for old in faces):raise ValueError('duplicate_topological_face')
            faces.append(f)
    N=np.array(rows);edges=tuple(sorted({tuple(sorted((a,b))) for f in faces for a,b in zip(f,f[1:]+f[:1])}));D2=np.zeros((len(edges),len(faces)),int);D1=np.zeros((len(vertices),len(edges)),int)
    for e,(a,b) in enumerate(edges):D1[a,e]=-1;D1[b,e]=1
    for j,f in enumerate(faces):
        for a,b in zip(f,f[1:]+f[:1]):e=edges.index(tuple(sorted((a,b))));D2[e,j]+=1 if (a,b)==edges[e] else -1
    if np.any(D1@D2) or np.any(np.count_nonzero(D2,axis=1)>2):raise ValueError('invalid_global_boundary')
    if any(abs(row.sum())>0 for row in D2[np.count_nonzero(D2,axis=1)==2]):raise ValueError('inconsistent_shared_orientation')
    T=[];seen=set()
    for a,b in te_edges:
        if a not in vertices or b not in vertices:raise ValueError('unknown_TE_vertex')
        aa,bb=vertices[a],vertices[b];key=tuple(sorted((aa,bb)))
        if key not in edges or key in seen:raise ValueError('invalid_TE_edge')
        seen.add(key);e=edges.index(key)
        if np.count_nonzero(D2[e])!=1:raise ValueError('internal_seam_is_not_TE')
        sign=1 if (aa,bb)==key else -1;T.append(-sign*D2[e])
    return dict(keys=keys,vertex_ids=tuple(ids),N=N,faces=faces,edges=edges,D1=D1,D2=D2,T=np.array(T))
