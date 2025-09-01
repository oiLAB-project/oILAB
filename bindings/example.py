import pyoilab as gb
import numpy as np
from scipy.spatial.transform import Rotation

# construct lattice
A= np.array([[0.5,0.5,0.0],
             [0.5,0.0,0.5],
             [0.0,0.5,0.5]])
lattice=gb.Lattice3D(A)

# construct a lattice vector from its cartesian coordinates
lv1 = lattice.latticeVector(np.array([1.0,2.0,2.0]))

# print its integer and cartesian coordinates
print("integer coordinates of lv1 = ", lv1.integerCoordinates())
cartesianCoordinates= lv1.cartesian()
print("cartesian coordinates of lv1 = ", cartesianCoordinates)

# reconstruct the lattice vector using the above cartesian coordinates
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

# Demonstrate the lattice box function
bv1= gb.LatticeVector3D(np.array([1,-1,0]),lattice)
bv2= gb.LatticeVector3D(np.array([1,1,-2]),lattice)
bv3= gb.LatticeVector3D(np.array([1,1,1]),lattice)
bv1= 3*bv1
bv2= 3*bv2
bv3= 3*bv3
config= lattice.box([bv1,bv2,bv3],"config.txt")

# Construct bicrystals with tilt axis 
misorientationAxis = lattice.reciprocalLatticeDirection(np.array([1.0,1.0,1.0]))
print("Rotation axis = ", misorientationAxis.cartesian())
rotations= lattice.generateCoincidentLattices(misorientationAxis)

for i, rotation in enumerate(rotations):
    rot = Rotation.from_matrix(rotation)
    angle = rot.magnitude() / 2
    axis = rot.as_rotvec()
    if np.allclose(axis, 0):
        angle=0.0
        axis=np.array([1,0,0])

    # Normalize axis and scale angle
    half_rotvec = axis / np.linalg.norm(axis) * angle
    positiveR = Rotation.from_rotvec(half_rotvec).as_matrix()
    negativeR = Rotation.from_rotvec(-half_rotvec).as_matrix()
    
    # From the two lattices by rotating the global lattice by +\theta/2 and -\theta/2
    lattice1 = gb.Lattice3D(A,positiveR)
    lattice2 = gb.Lattice3D(A,negativeR)

    try:
        bicrystal = gb.BiCrystal3D(lattice1,lattice2,False)
        print(f"Matrix {i}: angle - {np.degrees(2*angle):.2f} degrees, sigma = {bicrystal.sigma}")
        if (bicrystal.sigma<100):
            csl = bicrystal.csl
            # box vector along the axis
            bcslv1 = csl.latticeDirection(misorientationAxis.cartesian()).latticeVector()
            # box vector orthogonal to bcslv1
            #bcslv2 = csl.latticeDirection(np.array([1.0,0.0,-1.0])).latticeVector()
            bcslv2 = csl.latticeDirection(bcslv1.cross().cartesian()).latticeVector()
            # box vector orthogonal to bcslv1 and bcslv2
            bcslv3 = csl.latticeDirection(bcslv1.cross(bcslv2).cartesian()).latticeVector()
            config = bicrystal.box([bcslv1,bcslv2,bcslv3],1,1,"bc"+str(i))
    except Exception as e:
        print(f"Caught an exception: {e}")


