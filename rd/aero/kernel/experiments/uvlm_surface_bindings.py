"""Fixed material aero-point maps; PG is the structural implementation.

Q1 below is standalone AERO geometry interpolation, not a replacement shell.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import connection_shell as pg

def coordinates(m,n):
    xi=(1-np.cos(np.linspace(0,np.pi,m+1)))/2;eta=(1-np.cos(np.linspace(0,np.pi,n+1)))/2;xv=np.r_[xi[:-1]+np.diff(xi)/4,1.];xc=xi[:-1]+3*np.diff(xi)/4;yc=(eta[:-1]+eta[1:])/2
    return np.array([[x,y] for x in xv for y in eta]),np.array([[x,y] for x in xc for y in yc])

def q1(sid,m,n,corners):
    vertices,collocation=coordinates(m,n)
    def weights(p):
        x,y=p.T;N=np.column_stack(((1-x)*(1-y),(1-x)*y,x*y,x*(1-y)));Dx=np.column_stack((-(1-y),-y,y,1-y));Dy=np.column_stack((-(1-x),1-x,x,-x));return N,Dx,Dy
    Nv=weights(vertices)[0];Nc,Dx,Dy=weights(collocation);return dict(sid=sid,m=m,n=n,Nv=Nv,Nc=Nc,Dx=Dx,Dy=Dy,q0=np.array(corners,float))

def pg_shell(sid,m,n,model,q0):
    vertices,collocation=coordinates(m,n);Nv=pg.fields(model,np.column_stack((vertices,np.zeros(len(vertices)))))[0];p=np.column_stack((collocation,np.zeros(len(collocation))));Nc,G,_,_,X=pg.fields(model,p);J=pg.shell.geometry.reference(p[:,:2],p[:,2],model['h'],model['parameters'],False)[1];Dx=np.einsum('pib,pb->pi',G,J[:,:,0]);Dy=np.einsum('pib,pb->pi',G,J[:,:,1]);return dict(sid=sid,m=m,n=n,Nv=Nv,Nc=Nc,Dx=Dx,Dy=Dy,q0=np.asarray(q0))
