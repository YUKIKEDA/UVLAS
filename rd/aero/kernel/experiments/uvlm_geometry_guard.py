"""Whole-facet contact checks with allowed shared entities, no point sampling."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np

def cross2(a,b):return a[0]*b[1]-a[1]*b[0]
def inside(p,t,tol):
    a=t[1]-t[0];b=t[2]-t[0];q=p-t[0];det=cross2(a,b);u=cross2(q,b)/det;v=cross2(a,q)/det;return u>=-tol and v>=-tol and u+v<=1+tol
def segment_distance(p,a,b):
    e=b-a;s=np.dot(p-a,e)/np.dot(e,e);return np.linalg.norm(p-a-np.clip(s,0.,1.)*e)

def intersections(a,b,tol):
    na=np.cross(a[1]-a[0],a[2]-a[0]);nb=np.cross(b[1]-b[0],b[2]-b[0]);na/=np.linalg.norm(na);nb/=np.linalg.norm(nb);da=(b-a[0])@na;db=(a-b[0])@nb
    if min(da)>tol or max(da)<-tol or min(db)>tol or max(db)<-tol:return []
    out=[]
    if max(abs(da))<=tol and max(abs(db))<=tol:
        axis=np.argmax(abs(na));aa=np.delete(a,axis,axis=1);bb=np.delete(b,axis,axis=1)
        for p,q in zip(a,aa):
            if inside(q,bb,tol):out.append(p)
        for p,q in zip(b,bb):
            if inside(q,aa,tol):out.append(p)
        for i in range(3):
            for j in range(3):
                p=aa[i];r=aa[(i+1)%3]-p;q=bb[j];s=bb[(j+1)%3]-q;det=cross2(r,s)
                if abs(det)>tol:
                    u=cross2(q-p,s)/det;v=cross2(q-p,r)/det
                    if -tol<=u<=1+tol and -tol<=v<=1+tol:out.append(a[i]+u*(a[(i+1)%3]-a[i]))
    else:
        for x,y,n in [(a,b,nb),(b,a,na)]:
            axis=np.argmax(abs(n));yy=np.delete(y,axis,axis=1)
            for i in range(3):
                p=x[i];q=x[(i+1)%3];dp=np.dot(p-y[0],n);dq=np.dot(q-y[0],n)
                if abs(dp)<=tol and inside(np.delete(p,axis),yy,tol):out.append(p)
                if dp*dq<0:
                    r=p+(q-p)*(dp/(dp-dq))
                    if inside(np.delete(r,axis),yy,tol):out.append(r)
    return out

def validate_facets(facets):
    """facets: iterable (FaceId, three stable VertexIds, xyz[3,3])."""
    facets=list(facets);names=set();points={}
    for name,ids,x in facets:
        if name in names or len(set(ids))!=3 or not np.isfinite(x).all():raise ValueError('duplicate_or_invalid_facet')
        names.add(name);L=max(np.linalg.norm(x[i]-x[j]) for i in range(3) for j in range(i))
        if L<=0 or np.linalg.norm(np.cross(x[1]-x[0],x[2]-x[0]))<=1e-12*L*L:raise ValueError('degenerate_facet')
        for ident,p in zip(ids,x):
            if ident in points and np.linalg.norm(points[ident]-p)>1e-10*L:raise ValueError('shared_vertex_geometry_mismatch')
            points[ident]=p
    checked=0
    for i,(name,ids,a) in enumerate(facets):
        for other,jds,b in facets[i+1:]:
            L=max(np.linalg.norm(x-y) for x in a for y in b);tol=1e-10*L
            if np.any(a.max(axis=0)+tol<b.min(axis=0)) or np.any(b.max(axis=0)+tol<a.min(axis=0)):continue
            common=set(ids)&set(jds)
            if len(common)==3:raise ValueError('duplicate_geometric_facet')
            checked+=1;origin=a[0];aa=(a-origin)/L;bb=(b-origin)/L;hits=intersections(aa,bb,1e-10)
            for hit in hits:
                p=origin+L*hit
                if not common:raise ValueError('nonincident_facet_contact')
                if len(common)==1:
                    if np.linalg.norm(p-points[next(iter(common))])>10*tol:raise ValueError('intersection_beyond_shared_vertex')
                else:
                    first,second=tuple(common)
                    if segment_distance(p,points[first],points[second])>10*tol:raise ValueError('intersection_beyond_shared_edge')
    return dict(facets=len(facets),candidate_pairs_checked=checked)
