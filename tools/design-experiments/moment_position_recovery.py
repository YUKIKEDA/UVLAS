"""Reference-fixed C0 recovery of every native moment; structured diagnostic mesh.

This is a position map, not a replacement constitutive gradient or dynamics model.
All arrays are scalar maps to global moment coordinates, applied componentwise.
"""
import itertools
import numpy as np
import moment_enriched_ancf as f

def modal(p,ids):
    v=[np.polynomial.legendre.legvander(p[:,a],max(k[a] for k in ids)) for a in range(3)]
    return np.column_stack([np.prod([v[a][:,k[a]]*np.sqrt(2*k[a]+1) for a in range(3)],axis=0) for k in ids])

def facevalue(face,p):return modal(p,f.c.IDS)@face['Mcoeff']
def points(axis,sign,t):
    out=np.zeros((len(t),3),dtype=t.dtype);out[:,axis]=sign;out[:,[a for a in range(3) if a!=axis]]=t
    return out

class Recovery:
    def __init__(self,states,n,order=13,zorder=7,clamp=False):
        assert len(states)==n*n
        self.states=states;self.n=n;self.order=order;self.zorder=zorder;self.clamp=clamp
        self.facet={};self.cellfaces=[];offset=0;local=[]
        for cell,s in enumerate(states):
            j,i=divmod(cell,n);local.append((offset,offset+s['ni']));offset+=s['ni'];keys={}
            for face in s['faces'][2:]:
                a=face['axis'];sg=face['sign'];key=(a,i+(sg==1) if a==0 else i,j+(sg==1) if a==1 else j)
                if key not in self.facet:self.facet[key]=dict(face=face,offset=None,cell=cell)
                keys[a,sg]=key
            self.cellfaces.append(keys)
        for item in self.facet.values():item['offset']=offset;offset+=item['face']['n']
        self.size=offset;eye=np.eye(offset);self.native=[]
        for cell,s in enumerate(states):
            U=np.zeros((s['nscalar'],offset));a,b=local[cell];U[:s['ni']]=eye[a:b]
            for face in s['faces'][2:]:
                item=self.facet[self.cellfaces[cell][face['axis'],face['sign']]];p=face['points'].copy();p[:,face['axis']]=item['face']['sign']
                T=face['hat'].T@((face['w']/4)[:,None]*facevalue(item['face'],p));o=item['offset']
                U[face['offset']:face['offset']+face['n']]=T@eye[o:o+item['face']['n']]
            self.native.append(U)
        self.eye=eye;self.vertices={};self.edges={};self.facefields={};self.volume=[]
        # Only shared vertical traces contribute to the skeleton.
        for i,j,z in itertools.product(range(n+1),range(n+1),[-1,1]):
            values=[]
            for key,item in self.facet.items():
                a,x,y=key
                if a==0 and x==i and y<=j<=y+1:p=np.array([[item['face']['sign'],2*(j-y)-1,z]],float)
                elif a==1 and y==j and x<=i<=x+1:p=np.array([[2*(i-x)-1,item['face']['sign'],z]],float)
                else:continue
                values.append(self.raw(item,p)[0])
            self.vertices[i,j,z]=np.zeros(offset) if clamp and i==0 else np.mean(values,axis=0)
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

    def raw(self,item,p):
        o=item['offset'];return facevalue(item['face'],p)@self.eye[o:o+item['face']['n']]

    def correct_edge(self,base,left,right,degree,zero=False):
        if zero:return lambda t:np.zeros((len(t),self.size),dtype=t.dtype)
        end=base(np.array([-1.,1.]));A=np.polynomial.legendre.legvander(np.array([-1.,1.]),degree+1)[:,degree:degree+2]
        coeff=np.linalg.solve(A,np.array([left,right])-end)
        return lambda t:base(t)+np.polynomial.legendre.legvander(t,degree+1)[:,degree:degree+2]@coeff

    def edge(self,cell,direction,fixed,t):
        j,i=divmod(cell,self.n)
        if direction==2:key=(2,i+(fixed[0]==1),j+(fixed[1]==1),0)
        elif direction==0:key=(0,i,j+(fixed[1]==1),fixed[2])
        else:key=(1,i+(fixed[0]==1),j,fixed[2])
        return self.edges[key](t)

    def vertex(self,cell,signs):
        j,i=divmod(cell,self.n);return self.vertices[i+(signs[0]==1),j+(signs[1]==1),signs[2]]

    def face_boundary(self,cell,axis,sign,p):
        axes=[a for a in range(3) if a!=axis];out=np.zeros((len(p),self.size),dtype=p.dtype)
        for a in axes:
            b=next(x for x in axes if x!=a)
            for sg in [-1,1]:out+=((1+sg*p[:,b])/2)[:,None]*self.edge(cell,a,{axis:sign,b:sg},p[:,a])
        for x,y in itertools.product([-1,1],repeat=2):
            signs={axis:sign,axes[0]:x,axes[1]:y};blend=(1+x*p[:,axes[0]])*(1+y*p[:,axes[1]])/4
            out-=blend[:,None]*self.vertex(cell,signs)
        return out

    def make_face(self,cell,face):
        axis,sign=face['axis'],face['sign'];p,w=f.old.face_points(axis,sign,self.order,self.zorder);axes=[a for a in range(3) if a!=axis]
        M=facevalue(face,p);bubble=np.prod(1-p[:,axes]**2,axis=1);boundary=self.face_boundary(cell,axis,sign,p)
        U=self.native[cell][face['offset']:face['offset']+face['n']];gram=M.T@((w/4*bubble)[:,None]*M);coeff=np.linalg.solve(gram,U-M.T@((w/4)[:,None]*boundary))
        def evaluate(p):return self.face_boundary(cell,axis,sign,p)+np.prod(1-p[:,axes]**2,axis=1)[:,None]*(facevalue(face,p)@coeff)
        return evaluate

    def boundary(self,cell,p):
        out=np.zeros((len(p),self.size),dtype=p.dtype)
        for axis,sign in itertools.product(range(3),[-1,1]):
            key=('tb',cell,sign) if axis==2 else self.cellfaces[cell][axis,sign]
            # A shared face uses the representative's fixed coordinate.
            z=p.copy();z[:,axis]=sign if axis==2 else self.facet[key]['face']['sign']
            out+=((1+sign*p[:,axis])/2)[:,None]*self.facefields[key](z)
        for direction in range(3):
            axes=[a for a in range(3) if a!=direction]
            for x,y in itertools.product([-1,1],repeat=2):
                blend=(1+x*p[:,axes[0]])*(1+y*p[:,axes[1]])/4
                out-=blend[:,None]*self.edge(cell,direction,{axes[0]:x,axes[1]:y},p[:,direction])
        for signs in itertools.product([-1,1],repeat=3):out+=np.prod((1+p*np.array(signs))/2,axis=1)[:,None]*self.vertex(cell,dict(enumerate(signs)))
        return out

    def evaluate(self,cell,p):return self.boundary(cell,p)+np.prod(1-p*p,axis=1)[:,None]*(modal(p,f.VID)@self.volume[cell])

    def gradient(self,cell,p):
        vals=[]
        for axis in range(3):
            z=p.astype(complex);z[:,axis]+=1e-25j;vals.append(self.evaluate(cell,z).imag/1e-25)
        return np.stack(vals,axis=1)
