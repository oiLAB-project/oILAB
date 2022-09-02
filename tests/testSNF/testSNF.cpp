#include <LatticeModule.h>
#include <TextFileParser.h>

using namespace gbLAB;
int main()
{
    const auto A(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("A",true));
    const auto B(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("B",true));
    const auto F(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("F",true));

    Lattice<2> L1(A);
    Lattice<2> L2(B,F);
    BiCrystal<2> bc(L1,L2);
    std::cout<<"sigma_A="<<bc.sigmaA<<std::endl;
    std::cout<<"sigma_B="<<bc.sigmaB<<std::endl;
    if (std::fabs(bc.sigmaA) != 34 || std::fabs(bc.sigmaB) != 40)
    {
        throw std::runtime_error("SNF error in sigma calculation \n");
    }
    return 0;
}
