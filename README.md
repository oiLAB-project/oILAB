# oILAB

To install

1) mkdir build
2) cmake ..
3) make


Unit tests
----------

1) demonstrate lattice and reciprocal vectors, and directions
2) relationship between indices of a given vector w.r.t the four lattices in a bicrystal
3) get all lattices planes with period <= maxPeriod


src
----
1) Bicrystal class should contain the shift tensor \Lambda that gets computed when the bicrystal is constructed
2) introduce a new class Gb. A Gb is a bycrystal i.e.
            Class Gb : public Bicrystal
   It is constructed using a GB plane that inputted in terms of either the Cartesian coordinates, or its integer coordinates w.r.t to any one of the four lattices
3) Gb class contains a method that enumerates dicconnection modes for a given arbitrary DSC Burgers vector
