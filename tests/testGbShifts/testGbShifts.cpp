#include <TextFileParser.h>
#include <GbShifts.h>

// Hunter GB Details
//
// Sigma 3 [110]; theta = 70.52878
// GB normal: {1 - 1 2}
using namespace gbLAB;
int main()
{
    /*! [Types] */
    using VectorDimI = LatticeCore<3>::VectorDimI;
    using VectorDimD = LatticeCore<3>::VectorDimD;
    using Vector2d= Eigen::Vector2d;
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

    try
    {
        //double theta= 53.594515175286005615*std::numbers::pi/180;
        //double theta= 6.0089831977661480877*std::numbers::pi/180; // axis = 111
        double theta= 70.52878*std::numbers::pi/180;
        Eigen::AngleAxis<double> halfRotation(theta/2,rAxisGlobal.cartesian().normalized());
        Lattice<3> latticeA(lattice.latticeBasis,halfRotation.matrix());
        Lattice<3> latticeB(lattice.latticeBasis,halfRotation.matrix().transpose());

        BiCrystal<3> bc(latticeA,latticeB,false);

        // Specify GB normal
        VectorDimD gbNormal;
        //gbNormal << 17,23,-40;
        //gbNormal << 0.0, 0.0, 1.0;
        gbNormal << 1.0, 1.0, 2.0;
        gbNormal= halfRotation.matrix() * gbNormal;
        ReciprocalLatticeDirection<3> rd= latticeA.reciprocalLatticeDirection(gbNormal);


        // construct gb and shifts
        Gb<3> gb(bc,rd);
        ReciprocalLatticeVector<3> rAxisA(latticeA.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        LatticeVector<3> axisA(gb.bc.A.latticeDirection(axis).latticeVector());
        LatticeVector<3> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());
        auto basis= bc.csl.planeParallelLatticeBasis(
                bc.getReciprocalLatticeDirectionInC(gb.nA.reciprocalLatticeVector()), true);

        std::vector<LatticeVector<3>> gbCslVectors;
        gbCslVectors.push_back(gb.getPeriodVector(rAxisA));
        std::cout << "length of the period vector" << gb.getPeriodVector(rAxisA).cartesian().norm() << std::endl;
        gbCslVectors.push_back(axisC);
        GbShifts<3> shifts(gb,rAxisA,gbCslVectors,1.2);

        // construct the bicrystal
        std::vector<LatticeVector<3>> cslVectors;
        cslVectors.push_back(5*basis[0].latticeVector());
        cslVectors.push_back(4*gbCslVectors[0]);
        cslVectors.push_back(4*gbCslVectors[1]);
        auto points= gb.bc.box(cslVectors,0.5,1);
        gb.box(cslVectors,1,1,"ex.txt");


        std::ofstream config;

        int count= 0;
        for(const auto& pair : shifts.bShiftPairs)
        {
            std::string filename= "translate" + std::to_string(count) + ".txt";
            config.open(filename);
            config << points.size() << std::endl;
            config << "Lattice=\"";

            config << std::setprecision(15) << (2*cslVectors[0].cartesian()).transpose() << " ";
            config << std::setprecision(15) << (cslVectors[1].cartesian()).transpose() << " ";
            config << std::setprecision(15) << (cslVectors[2].cartesian()).transpose();
            config << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
            config << std::setprecision(15) << (-1 * cslVectors[0].cartesian()).transpose() << "\"" << std::endl;
            for (auto point : points)
            {
                if (&point.lattice == &bc.A)
                    config << "1 " << (point.cartesian() + pair.first.cartesian()/2).transpose() << " 0.05" << std::endl;
                if (&point.lattice == &bc.B)
                    config << "2 " << (point.cartesian() - pair.first.cartesian()/2).transpose() << " 0.05" << std::endl;
                if (&point.lattice == &bc.csl)
                    config << "3 " << (point.cartesian() + pair.second).transpose() << " 0.2" << std::endl;
                if (&point.lattice == &bc.dscl)
                    config << "4 " << point.cartesian().transpose() << " 0.01" << std::endl;
            }
            config.close();
            count++;
        }


    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Moving on the the next misorientation" << std::endl;
        std::cout << "-----------------------------------------------------------------------------"
                  << std::endl;
    }
    return 0;
}
