#include <LatticeModule.h>
#include <TextFileParser.h>

/* In this example, 
Lattice A = 
  0 0.5 0.5
0.5   0 0.5
0.5 0.5   0
maxis=
1
1
1
Cartesian coordinates of axis = 
1 1 1
###################################################
Misorientation angle = 6.0089831977661480877; Sigma = 273

Lattice B = 
0.001831501831501824995  0.52930402930402931094  0.46886446886446891957
 0.46886446886446891957 0.001831501831501824995  0.52930402930402931094
 0.52930402930402931094  0.46886446886446891957 0.001831501831501824995

Parallel CSL basis Cp= 
                    272.5                      -7.5                    -136.5
                      145                        -8                    -136.5
                    128.5                       0.5 8.8817841970012523234e-16

Parallel DSCL basis Dp = 
      0.998168498168498175                       -7.5                        0.5
    0.53113553113553113594                         -8                        0.5
    0.47069597069597068906                        0.5 -3.2534008047623634856e-18

Reduced DSCL basis vectors:
d1 =  0.001831501831501824995 -0.031135531135531135938  0.029304029304029310943
Integer coordinates of d1:-1  1 17

d2 =   0.02930402930402919992  0.001831501831501824995 -0.031135531135531024916
Integer coordinates of d2:-16  15 257

Reduced shift vectors: 
s1 = 0.0018315018315023801065    0.4688644688644695302   0.52930402930402986605
s2 =   0.52930402930402831174 0.0018315018315036013519   0.46886446886445387605

GBs of varying inclination (measured with respect to the first grain boundary)
-----------------------------------------------------------------------------
-- CAUTION: The integer coordinates of GB normals are w.r.t the reciprocal --
--          basis of the primitive unit cell.                              --
-----------------------------------------------------------------------------
1) Inclination = 0
nA =   5   6 -11
nB =   6   5 -11
GB period = 11.68332144554794283
CSL plane distance (Height)= 3.3726843908079144896
Glide disconnection Burgers vector =   0.031135531135532801272   -0.02930402930402919992 -0.0018315018315027131734; norm = 0.042796049251092468935
Step height of glide disconnection = 0.40768712416676056165
Step height of non-glide disconnection 1= 0.22237479499852241815
Step height of non-glide disconnection 2= 0.18531232916825057799
-----------------------------------------------------------------------------
2) Inclination = 1.9451190845039707522
nA =  17  23 -40
nB =  145  136 -281
GB period = 99.36548696604899078
CSL plane distance (Height)= 0.39655776945623399943
Glide disconnection Burgers vector =  0.80769230769237765344  -0.7307692307690558664 -0.07692307692309441336; norm = 1.0919284281982737372
Step height of glide disconnection = 0.091513331645593942731
Step height of non-glide disconnection 1= -0.18593184794961481465
Step height of non-glide disconnection 2= 0.19755258845548662183

*/


using namespace gbLAB;
int main()
{
    const int dim=3;
    /*! [lattice1] */
    Eigen::Matrix3d A;
    A << 0.0, 0.5, 0.5,
         0.5, 0.0, 0.5,
         0.5, 0.5, 0.0;
    Lattice<dim> L1(A);
    /*! [lattice1] */

    /*! [axis] */
    Eigen::Vector3d axis;
    axis << 1,1,1;
    ReciprocalLatticeVector<dim> rv(L1.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
    /*! [axis] */

    /*! [lattice2] */
    double angle= 6.0089831977661480877*M_PI/180;
    Eigen::AngleAxisd rotation(angle,axis.normalized());
    Lattice<dim> L2(A,rotation.matrix());
    /*! [lattice2] */

    try
    {
        /*! [bicrystal] */
        BiCrystal<dim> bc(L1,L2);
        std::cout << "Sigma = " << bc.sigma << std::endl;
        /*! [bicrystal] */

        /*! [gb] */
        ReciprocalLatticeVector<dim> normal(L1);
        normal << 17,23,-40;
        //normal << 5,6,-11;
        Gb<dim> gb(bc,normal);
        /*! [gb] */
        std::cout << "Miller indices w.r.t A: ";
        std::cout << gb.nA << std::endl;
        std::cout << "Miller indices w.r.t B: ";
        std::cout << gb.nB << std::endl;

        /*! [period vector] */
        LatticeVector<dim> glideA(rv.cross(gb.nA.reciprocalLatticeVector()).latticeVector());
        LatticeVector<dim> periodC(bc.getLatticeDirectionInC(glideA).latticeVector());
        /*! [period vector] */

        /*! [normal vector] */
        auto basis= bc.csl.planeParallelLatticeBasis(gb.bc.getReciprocalLatticeDirectionInC(gb.nA.reciprocalLatticeVector()),true);
        LatticeVector<dim> nonParallelC(basis[0].latticeVector());
        /*! [normal vector] */

        /*! [axis vector] */
        LatticeVector<dim> vectorAlongAxisA(bc.A.latticeVector(rv.cartesian()));
        LatticeVector<dim> vectorAlongAxisC(bc.getLatticeDirectionInC(vectorAlongAxisA).latticeVector());
        /*! [axis vector] */

        /*! [box vectors] */
        std::vector<LatticeVector<dim>> boxVectors;
        boxVectors.push_back(nonParallelC);
        boxVectors.push_back(periodC);
        boxVectors.push_back(vectorAlongAxisC);
        gb.box(boxVectors,0.8,2,"gb.txt",true);
        /*! [box vectors] */

    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        return -1;
    }
    return 0;
}
