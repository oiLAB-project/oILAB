"""
Python bindings for the oILAB C++ library
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['Lattice2D', 'Lattice3D', 'LatticeDirection2D', 'LatticeDirection3D', 'LatticeVector2D', 'LatticeVector3D', 'ReciprocalLatticeDirection2D', 'ReciprocalLatticeDirection3D', 'ReciprocalLatticeVector2D', 'ReciprocalLatticeVector3D']
class Lattice2D:
    def __init__(self, A: numpy.ndarray[numpy.float64[2, 2]], Q: numpy.ndarray[numpy.float64[2, 2]] = ...) -> None:
        ...
    def box(self, boxVectors: list[...], filename: str = '') -> list[...]:
        ...
    def generateCoincidentLattices(self, maxStrain: float, maxDen: float = 50, N: int = 30) -> list[numpy.ndarray[numpy.float64[2, 2]]]:
        ...
    def interPlanarSpacing(self, arg0: ...) -> float:
        ...
    def latticeVector(self, arg0: numpy.ndarray[numpy.float64[2, 1]]) -> ...:
        ...
    @property
    def F(self) -> numpy.ndarray[numpy.float64[2, 2]]:
        ...
    @property
    def latticeBasis(self) -> numpy.ndarray[numpy.float64[2, 2]]:
        ...
    @property
    def reciprocalBasis(self) -> numpy.ndarray[numpy.float64[2, 2]]:
        ...
class Lattice3D:
    def __init__(self, A: numpy.ndarray[numpy.float64[3, 3]], Q: numpy.ndarray[numpy.float64[3, 3]] = ...) -> None:
        ...
    def box(self, boxVectors: list[...], filename: str = '') -> list[...]:
        ...
    def generateCoincidentLattices(self, rd: ..., maxDen: float = 100, N: int = 100) -> list[numpy.ndarray[numpy.float64[3, 3]]]:
        ...
    def interPlanarSpacing(self, arg0: ...) -> float:
        ...
    def latticeVector(self, arg0: numpy.ndarray[numpy.float64[3, 1]]) -> ...:
        ...
    @property
    def F(self) -> numpy.ndarray[numpy.float64[3, 3]]:
        ...
    @property
    def latticeBasis(self) -> numpy.ndarray[numpy.float64[3, 3]]:
        ...
    @property
    def reciprocalBasis(self) -> numpy.ndarray[numpy.float64[3, 3]]:
        ...
class LatticeDirection2D:
    @typing.overload
    def __init__(self, arg0: LatticeVector2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: LatticeDirection2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.int64[2, 1]], arg1: Lattice2D) -> None:
        ...
    def cartesian(self) -> numpy.ndarray[numpy.float64[2, 1]]:
        """
        Cartesian coordinates of the lattice direction.
        """
    def dot(self, arg0: ...) -> int:
        """
        dot product with a reciprocal lattice vector
        """
    def integerCoordinates(self) -> numpy.ndarray[numpy.int64[2, 1]]:
        """
        Integer coordinates of the lattice direction.
        """
    def latticeVector(self) -> LatticeVector2D:
        """
        Get the lattice vector
        """
class LatticeDirection3D:
    @typing.overload
    def __init__(self, arg0: LatticeVector3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: LatticeDirection3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.int64[3, 1]], arg1: Lattice3D) -> None:
        ...
    def cartesian(self) -> numpy.ndarray[numpy.float64[3, 1]]:
        """
        Cartesian coordinates of the lattice direction.
        """
    def dot(self, arg0: ...) -> int:
        """
        dot product with a reciprocal lattice vector
        """
    def integerCoordinates(self) -> numpy.ndarray[numpy.int64[3, 1]]:
        """
        Integer coordinates of the lattice direction.
        """
    def latticeVector(self) -> LatticeVector3D:
        """
        Get the lattice vector
        """
class LatticeVector2D:
    def __add__(self, arg0: LatticeVector2D) -> LatticeVector2D:
        ...
    def __iadd__(self, arg0: LatticeVector2D) -> LatticeVector2D:
        ...
    @typing.overload
    def __init__(self, arg0: Lattice2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.float64[2, 1]], arg1: Lattice2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.int64[2, 1]], arg1: Lattice2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: LatticeVector2D) -> None:
        ...
    def __isub__(self, arg0: LatticeVector2D) -> LatticeVector2D:
        ...
    def __mul__(self, arg0: int) -> LatticeVector2D:
        ...
    def __rmul__(self, arg0: int) -> LatticeVector2D:
        ...
    def __sub__(self, arg0: LatticeVector2D) -> LatticeVector2D:
        ...
    def cartesian(self) -> numpy.ndarray[numpy.float64[2, 1]]:
        ...
    @typing.overload
    def cross(self, arg0: LatticeVector2D) -> ...:
        ...
    @typing.overload
    def cross(self) -> ...:
        ...
    def dot(self, arg0: ...) -> int:
        """
        dot product with a reciprocal lattice vector
        """
    @typing.overload
    def integerCoordinates(self) -> numpy.ndarray[numpy.int64[2, 1]]:
        """
        output the integer coordinates of the lattice vecctor
        """
    @typing.overload
    def integerCoordinates(self, arg0: numpy.ndarray[numpy.int64[2, 1]]) -> None:
        """
        input the integer coordinates of the lattice vector
        """
class LatticeVector3D:
    def __add__(self, arg0: LatticeVector3D) -> LatticeVector3D:
        ...
    def __iadd__(self, arg0: LatticeVector3D) -> LatticeVector3D:
        ...
    @typing.overload
    def __init__(self, arg0: Lattice3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: Lattice3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.int64[3, 1]], arg1: Lattice3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: LatticeVector3D) -> None:
        ...
    def __isub__(self, arg0: LatticeVector3D) -> LatticeVector3D:
        ...
    def __mul__(self, arg0: int) -> LatticeVector3D:
        ...
    def __rmul__(self, arg0: int) -> LatticeVector3D:
        ...
    def __sub__(self, arg0: LatticeVector3D) -> LatticeVector3D:
        ...
    def cartesian(self) -> numpy.ndarray[numpy.float64[3, 1]]:
        ...
    @typing.overload
    def cross(self, arg0: LatticeVector3D) -> ...:
        ...
    @typing.overload
    def cross(self) -> ...:
        ...
    def dot(self, arg0: ...) -> int:
        """
        dot product with a reciprocal lattice vector
        """
    @typing.overload
    def integerCoordinates(self) -> numpy.ndarray[numpy.int64[3, 1]]:
        """
        output the integer coordinates of the lattice vecctor
        """
    @typing.overload
    def integerCoordinates(self, arg0: numpy.ndarray[numpy.int64[3, 1]]) -> None:
        """
        input the integer coordinates of the lattice vector
        """
class ReciprocalLatticeDirection2D:
    @typing.overload
    def __init__(self, arg0: ReciprocalLatticeVector2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: ReciprocalLatticeDirection2D) -> None:
        ...
    def cartesian(self) -> numpy.ndarray[numpy.float64[2, 1]]:
        """
        Cartesian coordinates of the reciprocal lattice direction.
        """
    def dot(self, arg0: LatticeVector2D) -> int:
        """
        dot product with a lattice vector
        """
    def integerCoordinates(self) -> numpy.ndarray[numpy.int64[2, 1]]:
        """
        Integer coordinates of the reciprocal lattice direction.
        """
    def planeSpacing(self) -> float:
        """
        distance between two consecutive planes
        """
    def reciprocalLatticeVector(self) -> ReciprocalLatticeVector2D:
        """
        Get the reciprocal lattice vector
        """
    def stacking(self) -> int:
        """
        stacking
        """
class ReciprocalLatticeDirection3D:
    @typing.overload
    def __init__(self, arg0: ReciprocalLatticeVector3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: ReciprocalLatticeDirection3D) -> None:
        ...
    def cartesian(self) -> numpy.ndarray[numpy.float64[3, 1]]:
        """
        Cartesian coordinates of the reciprocal lattice direction.
        """
    def dot(self, arg0: LatticeVector3D) -> int:
        """
        dot product with a lattice vector
        """
    def integerCoordinates(self) -> numpy.ndarray[numpy.int64[3, 1]]:
        """
        Integer coordinates of the reciprocal lattice direction.
        """
    def planeSpacing(self) -> float:
        """
        distance between two consecutive planes
        """
    def reciprocalLatticeVector(self) -> ReciprocalLatticeVector3D:
        """
        Get the reciprocal lattice vector
        """
    def stacking(self) -> int:
        """
        stacking
        """
class ReciprocalLatticeVector2D:
    def __add__(self, arg0: ReciprocalLatticeVector2D) -> ReciprocalLatticeVector2D:
        ...
    def __iadd__(self, arg0: ReciprocalLatticeVector2D) -> ReciprocalLatticeVector2D:
        ...
    @typing.overload
    def __init__(self, arg0: Lattice2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.float64[2, 1]], arg1: Lattice2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.int64[2, 1]], arg1: Lattice2D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: ReciprocalLatticeVector2D) -> None:
        ...
    def __isub__(self, arg0: ReciprocalLatticeVector2D) -> ReciprocalLatticeVector2D:
        ...
    def __mul__(self, arg0: int) -> ReciprocalLatticeVector2D:
        ...
    def __rmul__(self, arg0: int) -> ReciprocalLatticeVector2D:
        ...
    def __sub__(self, arg0: ReciprocalLatticeVector2D) -> ReciprocalLatticeVector2D:
        ...
    def cartesian(self) -> numpy.ndarray[numpy.float64[2, 1]]:
        ...
    def cross(self, arg0: ReciprocalLatticeVector2D) -> LatticeDirection2D:
        ...
    def dot(self, arg0: LatticeVector2D) -> int:
        ...
    def integerCoordinates(self) -> numpy.ndarray[numpy.int64[2, 1]]:
        ...
class ReciprocalLatticeVector3D:
    def __add__(self, arg0: ReciprocalLatticeVector3D) -> ReciprocalLatticeVector3D:
        ...
    def __iadd__(self, arg0: ReciprocalLatticeVector3D) -> ReciprocalLatticeVector3D:
        ...
    @typing.overload
    def __init__(self, arg0: Lattice3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: Lattice3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[numpy.int64[3, 1]], arg1: Lattice3D) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: ReciprocalLatticeVector3D) -> None:
        ...
    def __isub__(self, arg0: ReciprocalLatticeVector3D) -> ReciprocalLatticeVector3D:
        ...
    def __mul__(self, arg0: int) -> ReciprocalLatticeVector3D:
        ...
    def __rmul__(self, arg0: int) -> ReciprocalLatticeVector3D:
        ...
    def __sub__(self, arg0: ReciprocalLatticeVector3D) -> ReciprocalLatticeVector3D:
        ...
    def cartesian(self) -> numpy.ndarray[numpy.float64[3, 1]]:
        ...
    def cross(self, arg0: ReciprocalLatticeVector3D) -> LatticeDirection3D:
        ...
    def dot(self, arg0: LatticeVector3D) -> int:
        ...
    def integerCoordinates(self) -> numpy.ndarray[numpy.int64[3, 1]]:
        ...
