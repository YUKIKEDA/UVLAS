"""T-E24: tangent two-point edge traces times quadratic thickness."""
import itertools
import numpy as np
from q9_trace_shear import cardinal
A=1/np.sqrt(3.)
AXES=[([-A,A],[-1.,1.]),([-1.,1.],[-A,A])]
def interpolation(points):
    points=np.asarray(points);S=np.zeros((len(points),2,24));column=0
    for component,(xs,ys) in enumerate(AXES):
        for i,j,k in itertools.product(range(2),range(2),range(3)):
            S[:,component,column]=cardinal(points[:,0],xs,i)*cardinal(points[:,1],ys,j)*cardinal(points[:,2],[-1.,0.,1.],k);column+=1
    return S
def sampling_points():
    points=[];components=[]
    for component,(xs,ys) in enumerate(AXES):
        for x,y,z in itertools.product(xs,ys,[-1.,0.,1.]):points.append([x,y,z]);components.append(component)
    return np.array(points),np.array(components)
