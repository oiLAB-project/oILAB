#include <LatticeModule.h>
#include <TextFileParser.h>
#include <BiCrystal.h>

using namespace gbLAB;

int main()
{
    /*! [Types] */
    using VectorDimI = LatticeCore<3>::VectorDimI;
    using VectorDimD = LatticeCore<3>::VectorDimD;
    using IntScalarType = LatticeCore<3>::IntScalarType;
    /*! [Types] */

    /*! [Lattice] */
    const auto A(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("A",true));
    Lattice<3> lattice(A);
    std::cout << "Lattice A = " << std::endl;
    std::cout << lattice.latticeBasis << std::endl;
    /*! [Lattice] */

    /*! [Axis] */
    const auto axis (TextFileParser("bicrystal_3d.txt").readMatrix<double,3,1>("axis",true));
    ReciprocalLatticeVector<3> rAxisGlobal(lattice.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
    std::cout << "Cartesian coordinates of axis = " << std::endl;
    std::cout << rAxisGlobal.cartesian().transpose() << std::endl;
    /*! [Axis] */

    /*! [Generate bicrystal] */
    const auto& coincidentLattices= lattice.generateCoincidentLattices(rAxisGlobal,150,150);
    /*! [Generate bicrystal] */

    fstream surfaceEnergyFile;
    surfaceEnergyFile.open("surfaceEnergy.txt",ios::out);
    /*! [Misorientation] */
    for (const auto& rotation : coincidentLattices)
    {
        // Loop over misorientation angles
        /*! [Misorientation] */
        try
        {
            /*! [SNF] */
            double theta= acos((rotation.trace()-1.0)/2.0);
            Eigen::AngleAxis<double> halfRotation(theta/2,rAxisGlobal.cartesian().normalized());
            Lattice<3> latticeA(lattice.latticeBasis,halfRotation.matrix());
            Lattice<3> latticeB(lattice.latticeBasis,halfRotation.matrix().transpose());

            BiCrystal<3> bc(latticeA,latticeB,false);
            //if (bc.sigma > 20) continue;
            std::cout << "Misorientation angle = " << std::setprecision(20) << theta*180/M_PI << "; ";
            std::cout << "Sigma = " << std::setprecision(20) << bc.sigma << std::endl;
            std::cout << std::endl;
            std::cout << "Lattice A = " << std::endl;
            std::cout << std::setprecision(20) << latticeA.latticeBasis << std::endl;
            std::cout << "Lattice B = " << std::endl;
            std::cout << std::setprecision(20) << latticeB.latticeBasis << std::endl;
            std::cout << std::endl;
            std::cout << "Parallel CSL basis Cp= " << std::endl;
            std::cout << std::setprecision(20) << bc.csl.latticeBasis <<  std::endl;
            std::cout << std::endl;
            std::cout << "Parallel DSCL basis Dp = " << std::endl;
            std::cout << std::setprecision(20) << bc.dscl.latticeBasis <<  std::endl;
            std::cout << std::endl;
            /*! [SNF] */

            /*! [Invariance] */
            auto reducedDsclBasis= RLLL(bc.dscl.latticeBasis,0.75);
            auto U_Dscl= reducedDsclBasis.unimodularMatrix();

            std::cout << "Reduced DSCL basis vectors:" << std::endl;
            std::cout << "d1 = ";
            std::cout << std::setprecision(20) << reducedDsclBasis.reducedBasis().col(0).transpose() << std::endl;
            std::cout << "Integer coordinates of d1:";
            LatticeVector<3> d1(bc.dscl);
            d1 << U_Dscl.col(0).template cast<IntScalarType>();
            std::cout << std::setprecision(20) << d1.transpose() << std::endl;
            std::cout << std::endl;

            LatticeVector<3> d2(bc.dscl);
            std::cout << "d2 = ";
            std::cout << std::setprecision(20) << reducedDsclBasis.reducedBasis().col(1).transpose() << std::endl;
            std::cout << "Integer coordinates of d2:";
            d2 << U_Dscl.col(1).template cast<IntScalarType>();
            std::cout << std::setprecision(20) << d2.transpose() << std::endl;
            std::cout << std::endl;

            // Specify GB normal
            VectorDimD gbNormal;
            gbNormal << 0.0, 0.0, 1.0;
            ReciprocalLatticeDirection<3> rd= latticeA.reciprocalLatticeDirection(gbNormal);

            Gb<3> gb(bc,rd);
            ReciprocalLatticeVector<3> rAxisA(latticeA.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
            try
            {
                LatticeVector<3> glideA(rAxisA.cross(gb.nA.reciprocalLatticeVector()).latticeVector());
                LatticeVector<3> burgersVectorA(bc.getLatticeDirectionInD(glideA).latticeVector());
                LatticeVector<3> periodVectorA(bc.getLatticeDirectionInC(glideA).latticeVector());

                ReciprocalLatticeVector<3> gbNormalA = gb.nA.reciprocalLatticeVector();
                ReciprocalLatticeVector<3> gbNormalCsl= bc.getReciprocalLatticeDirectionInC(
                        gbNormalA).reciprocalLatticeVector();

                /*! [Output] */
                double epsilon = 1e-8;
                if (periodVectorA.cartesian().norm() < 50) {
                    std::cout << "nA = " << gb.nA << std::endl;
                    std::cout << "nB = " << gb.nB << std::endl;
                    std::cout << "GB period = " << std::setprecision(20) << periodVectorA.cartesian().norm()
                              << std::endl;
                    std::cout << "CSL plane distance (Height)= " << std::setprecision(20)
                              << 1.0 / gbNormalCsl.cartesian().norm() << std::endl;

                    // planes perpendicular to the GB
                    Eigen::Vector3d doubleCoordinates= latticeA.latticeBasis.transpose()*burgersVectorA.cartesian();
                    auto integerCoordinates= LatticeCore<3>::rationalApproximation(doubleCoordinates);
                    ReciprocalLatticeVector<3> temp(integerCoordinates,latticeA);
                    ReciprocalLatticeDirection<3> planeOrthogonalToGb(temp);
                    int numberOfPlanesOrthogonalToGB= planeOrthogonalToGb.stacking();
                    std::cout << "Number of lattice planes perpendicular to the GB = "
                              << numberOfPlanesOrthogonalToGB << std::endl;
                    std::cout << "Number of lattice planes parallel to the GB = "
                              << gb.nA.stacking() << std::endl;



                    LatticeVector<3> axisCSL(bc.csl.latticeDirection(axis).latticeVector());
                    auto basis= latticeA.planeParallelLatticeBasis(gb.nA);
                    int numberOfCslPoints= std::round(
                            axisCSL.cartesian().cross(periodVectorA.cartesian()).norm()/
                            basis[1].cartesian().cross(basis[2].cartesian()).norm());
                    std::cout << "Number of CSL points = " << numberOfCslPoints << std::endl;
                    surfaceEnergyFile << theta*180/M_PI  << "\t"
                                      << bc.sigma << "\t"
                                      << numberOfCslPoints << "\t"
                                      << axisCSL.cartesian().cross(periodVectorA.cartesian()).norm()  << "\t"
                                      << numberOfCslPoints/axisCSL.cartesian().cross(periodVectorA.cartesian()).norm()
                                      << std::endl;

                    IntScalarType alpha= basis[1].dot(planeOrthogonalToGb);
                    IntScalarType beta= basis[2].dot(planeOrthogonalToGb);
                    IntScalarType x,y;
                    IntegerMath<IntScalarType>::extended_gcd(alpha, beta, x, y);
                    LatticeVector<3> gbBasisAtom(x * basis[1].latticeVector() + y * basis[2].latticeVector());
                    std::vector<VectorDimI> basisAtomIndex;

                    for(int j=0; j< gb.nA.stacking(); ++j) {
                        for (int i = 0; i < numberOfCslPoints; ++i) {
                            LatticeVector<3> atom((i+1)*gbBasisAtom);
                            atom= atom + j*basis[0].latticeVector();
                            VectorDimI index;
                            index(0) =
                                IntegerMath<IntScalarType>::positive_modulo(atom.dot(planeOrthogonalToGb) , numberOfPlanesOrthogonalToGB);
                            index(1) =
                                IntegerMath<IntScalarType>::positive_modulo(atom.dot(rAxisA) , (ReciprocalLatticeDirection<3>(rAxisA).stacking()));
                            index(2) = j;
                            basisAtomIndex.push_back(index);
                            //std::cout << index.transpose() << std::endl;
                        }
                   }


                    /*
                    std::vector<VectorDimI> basisAtomsIntegerCoordinates;
                    int count= 0;
                    for(int i=0;i<numberOfCslPoints; ++i) {
                        for (int j = 0; j < numberOfCslPoints; ++j) {
                            LatticeVector<3> v(latticeA);
                            v= i*basis[1].latticeVector() + j*basis[2].latticeVector();
                            if(v.dot(planeOrthogonalToGb) < 0) v=-1*v;
                            v = v - (v.dot(planeOrthogonalToGb) / numberOfPlanesOrthogonalToGB) *
                                    LatticeVector<3>(periodVectorA.cartesian(), latticeA);
                            basisAtomsIntegerCoordinates.push_back((VectorDimI) v);
                            count ++;
                            if (count == numberOfCslPoints) break;
                        }
                    }
                    if (count != numberOfCslPoints) exit(0);

                    for(auto vec : basisAtomsIntegerCoordinates)
                        std::cout << vec.dot(planeOrthogonalToGb.reciprocalLatticeVector()) << std::endl;
                        */


                    std::map<int,int> mesostate;
                    for(auto index : basisAtomIndex)
                    {
                        mesostate.insert(std::make_pair<int,int>(index(0),index(2)));
                    }
                    for(auto pair : mesostate)
                        std::cout << pair.first << "  " << pair.second << std::endl;

                    // mesostate generator
                    std::vector<std::map<int,int>> mesostates;
                    mesostates.push_back(mesostate);


                    std::cout << "-----------------------------------------------------------------------------"
                              << std::endl;
                }
            }
                /*! [Output] */
            catch(std::runtime_error& e)
            {
                std::cout << e.what() << std::endl;
                std::cout << "Moving onto the next inclination" << std::endl;
            }
        }
        catch(std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
            std::cout << "Moving on the the next misorientation" << std::endl;
            std::cout << "-----------------------------------------------------------------------------"
                      << std::endl;
        }
    }

    surfaceEnergyFile.close();
    return 0;
}
