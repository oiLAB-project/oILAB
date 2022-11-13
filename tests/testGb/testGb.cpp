#include <LatticeModule.h>
#include <TextFileParser.h>

using namespace gbLAB;
int main()
{
    /*! [Lattice] */
    const auto A(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("A",true));
    const auto B(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("B",true));
    const auto F(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("F",true));

    Lattice<2> L1(A);
    Lattice<2> L2(B,F);
    /*! [Lattice] */

    /*! [GB] */
    try
    {
        BiCrystal<2> bc(L1, L2);
        ReciprocalLatticeVector<2> rvec(L2);
        rvec << 1,1;
        ReciprocalLatticeDirection<2> n(rvec);
        Gb<2> gb(bc,n);
        std::cout << "Miller indices w.r.t A: ";
        std::cout << gb.nA << std::endl;
        std::cout << "Miller indices w.r.t B: ";
        std::cout << gb.nB << std::endl;

        LatticeVector<2> b(bc.dscl);
        b << 1,1;
        std::cout << "Step height A: ";
        double stepHeightA= gb.stepHeightA(b);
        std::cout << stepHeightA << std::endl;
        std::cout << "Step height B: ";
        double stepHeightB= gb.stepHeightB(b);
        std::cout << stepHeightB << std::endl;

        if((gb.nA.cartesian().normalized()+gb.nB.cartesian().normalized()).norm() > 1e-6)
            throw std::runtime_error("Cartesian coordinates of normalized nA and nB are not anti-parallel");
        if(abs(stepHeightA-stepHeightB - b.cartesian().dot(gb.nA.cartesian().normalized())) > 1e-6)
            throw std::runtime_error("Error in step height calculation");
    }
    /*! [GB] */
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        return -1;
    }
    return 0;
}
