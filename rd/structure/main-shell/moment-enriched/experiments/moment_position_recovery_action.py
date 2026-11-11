"""Apply the same C0 recovery to moment values, without a global identity."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
import moment_position_recovery as old
f=old.f;modal=old.modal;facevalue=old.facevalue;points=old.points
def light(s):
    keys=['axis','sign','offset','n','Mcoeff','completion','points','w','hat']
    return dict(nscalar=s['nscalar'],ni=s['ni'],faces=[{key:face[key] for key in keys} for face in s['faces']])
class RecoveryAction(old.Recovery):
    def __init__(self,states,n,values,order=13,zorder=7,clamp=False):
        assert len(states)==n*n and len(values)==n*n
        self.states=states;self.n=n;self.order=order;self.zorder=zorder;self.clamp=clamp
        self.native=[np.asarray(v) for v in values];self.size=self.native[0].shape[1]
        self.facet={};self.cellfaces=[];self.vertices={};self.edges={};self.facefields={};self.volume=[]
        for cell,(s,value) in enumerate(zip(states,self.native)):
            assert value.shape==(s['nscalar'],self.size)
            j,i=divmod(cell,n);keys={}
            for face in s['faces'][2:]:
                a,sg=face['axis'],face['sign'];key=(a,i+(sg==1) if a==0 else i,j+(sg==1) if a==1 else j);u=value[face['offset']:face['offset']+face['n']]
                if key not in self.facet:self.facet[key]=dict(face=face,values=u.copy(),cell=cell)
                else:
                    previous=self.facet[key];T=face['Mcoeff'].T@previous['face']['Mcoeff'];difference=np.linalg.norm(u-T@previous['values'])
                    if difference>1e-9*max(1.,np.linalg.norm(u)):raise ValueError('recovery_action_incompatible_shared_moments')
                keys[a,sg]=key
            self.cellfaces.append(keys)
        # Only shared vertical traces contribute to the skeleton.
        for i,j,z in itertools.product(range(n+1),range(n+1),[-1,1]):
            values=[]
            for key,item in self.facet.items():
                a,x,y=key
                if a==0 and x==i and y<=j<=y+1:p=np.array([[item['face']['sign'],2*(j-y)-1,z]],float)
                elif a==1 and y==j and x<=i<=x+1:p=np.array([[2*(i-x)-1,item['face']['sign'],z]],float)
                else:continue
                values.append(self.raw(item,p)[0])
            self.vertices[i,j,z]=np.zeros(self.size) if clamp and i==0 else np.mean(values,axis=0)
        # Vertical edge polynomials, with endpoint correction orthogonal to P4.
        for i,j in itertools.product(range(n+1),repeat=2):
            incident=[]
            for key,item in self.facet.items():
                a,x,y=key
                if a==0 and x==i and y<=j<=y+1:incident.append((item,2*(j-y)-1))
                elif a==1 and y==j and x<=i<=x+1:incident.append((item,2*(i-x)-1))
            def base(t,incident=incident):
                vals=[]
                for item,u in incident:
                    a=item['face']['axis'];p=np.zeros((len(t),3),dtype=t.dtype);p[:,a]=item['face']['sign'];p[:,1-a]=u;p[:,2]=t;vals.append(self.raw(item,p))
                return sum(vals)/len(vals)
            degree=5 if any(max(item['face']['completion']['ranks'])>0 for item,u in incident) else 3
            self.edges[2,i,j,0]=self.correct_edge(base,self.vertices[i,j,-1],self.vertices[i,j,1],degree,clamp and i==0)
        # Horizontal edges come from their single shared vertical trace.
        for key,item in self.facet.items():
            a,i,j=key;direction=1-a
            for z in [-1,1]:
                def base(t,item=item,z=z,direction=direction):
                    p=np.zeros((len(t),3),dtype=t.dtype);p[:,item['face']['axis']]=item['face']['sign'];p[:,2]=z;p[:,direction]=t;return self.raw(item,p)
                end=(i+(direction==0),j+(direction==1),z)
                degree=9 if max(item['face']['completion']['ranks'])>0 else 5
                self.edges[direction,i,j,z]=self.correct_edge(base,self.vertices[i,j,z],self.vertices[end],degree,clamp and i==0 and direction==1)
        # Build a shared vertical face once. Each cell owns its two TB faces.
        for cell,s in enumerate(states):
            j,i=divmod(cell,n)
            for face in s['faces']:
                a,sg=face['axis'],face['sign'];key=('tb',cell,sg) if a==2 else self.cellfaces[cell][a,sg]
                if key not in self.facefields:
                    self.facefields[key]=self.make_face(cell,face)
            p,w=f.volume_points(order,zorder);V=modal(p,f.VID);bubble=np.prod(1-p*p,axis=1);boundary=self.boundary(cell,p)
            gram=V.T@((w/8*bubble)[:,None]*V);rhs=self.native[cell][:f.NV]-V.T@((w/8)[:,None]*boundary)
            self.volume.append(np.linalg.solve(gram,rhs))

    def raw(self,item,p):return facevalue(item['face'],p)@item['values']
