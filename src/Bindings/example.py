import pyoilab as gb
import numpy as np

from ovito.io import import_file
from ovito.vis import Viewport
from ovito.data import Particles

# construct lattice
lattice=gb.Lattice3D(np.array([[0.5,0.5,0.0],
                               [0.5,0.0,0.5],
                               [0.0,0.5,0.5]]))

# construct a lattice vector from its cartesian coordinates
lv1 = lattice.latticeVector(np.array([1.0,2.0,2.0]))

# print its integer and cartesian coordinates
print("integer coordinates of lv1 = ", lv1.integerCoordinates())
cartesianCoordinates= lv1.cartesian()
print("cartesian coordinates of lv1 = ", cartesianCoordinates)

# reconstruct the lattice vectors using the above cartesian coordinates
lv2 = gb.LatticeVector3D(lv1.cartesian(),lattice)

# confirm that lv1 and lv2 are the same
#$print("integer coorindates of lv2 = ", lv2.integerCoordinates())
print("Are lv1 and lv2 identical: ", lv1.integerCoordinates()==lv2.integerCoordinates())

# now change lv2's integer coordinates
lv2.integerCoordinates(np.array([1,2,3],dtype=np.int64))
print("Integer coordinates of lv2  changed to = ", lv2.integerCoordinates())
print("Cartesian coordinates of lv2 = ", lv2.cartesian())

# Demonstrate addition and multiplication
lv3=4*(lv1+lv2)
print("Integer coordinates of lv3 = ", lv3.integerCoordinates())
print("Cartesian coordinates of lv3 = ", lv3.cartesian())


# Construct reciprocal lattice vectors
rlv1= gb.ReciprocalLatticeVector3D(np.array([1,5,6],dtype=np.int64),lattice)
rlv2= gb.ReciprocalLatticeVector3D(np.array([2.,2.,4.],dtype=np.float64),lattice)

# Demonstrate cross product
ld=rlv1.cross(rlv2)
print("Integer coorindates of the resulting lattice direction = ", ld.integerCoordinates())
print("Is ld perpendicular to rlv1 and rlv2?", (ld.dot(rlv1)==0 and ld.dot(rlv2)==0))

# Demonstrate box function
bv1= gb.LatticeVector3D(np.array([1,-1,0]),lattice)
bv2= gb.LatticeVector3D(np.array([1,1,-2]),lattice)
bv3= gb.LatticeVector3D(np.array([1,1,1]),lattice)
config= lattice.box([bv1,bv2,bv3],"config.txt")

rotations= lattice.generateCoincidentLattices(gb.ReciprocalLatticeDirection3D(rlv1))

def create(source,frame):
    pipeline = import_file("/Users/Nikhil/Documents/Academic/Software/oILAB/build/src/config.txt")
    pipeline.add_to_scene()
    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    return pipeline

