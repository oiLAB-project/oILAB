#include <LatticeModule.h>
#include <TextFileParser.h>
#include <range.h>

using namespace gbLAB;
int main()
{
    const int m_max=19; 
    const int n_max=19; 
    const auto A(TextFileParser("cubic_lattice.txt").readMatrix<double,3,3>("A",true));
    Lattice<3> L1(A);

    //LatticeDirection<3> axis(TextFileParser("cubic_lattice.txt").readMatrix<long long,3,1>("axis",true),L1);
    const auto axis(L1.latticeDirection(TextFileParser("cubic_lattice.txt").readMatrix<double,3,1>("axis",true)));
    //const auto unit_axis(axis.cartesian().normalized());
    const auto unit_axis(axis.cartesian().normalized());
    std::cout << unit_axis << std::endl;

    Eigen::Matrix<double,3,3> R;
    for (auto m : range<int>(1,m_max))
        for (auto n : range<int>(0,n_max))
        {
            m = m/IntegerMath<int>::gcd(m,n);
            n = n/IntegerMath<int>::gcd(m,n);
            // rotation that guarantees coincidence 
            // RANGANATHAN, S. (1966). Acta Cryst. A21, 197-199.
            double theta = 2.0*atan(axis.template cast<double>().norm()*n/m);
            int S = pow(m,2)+axis.squaredNorm() * pow(n,2);
            std::cout<<"(m,n)=("<<m<<","<<n<<"): "
                     <<", S="<<S
                     <<", theta="<<theta*180/M_PI << std::flush;



            R = Eigen::AngleAxis<double>(theta,unit_axis);


            // create bicrystal 
            Lattice<3> L2(A,R);
            BiCrystal<3> bc(L1,L2,false);

            std::cout<<", sigma="<<bc.sigma
                     <<", alpha="<<S/bc.sigma
                     <<std::endl; 
            if(S/bc.sigma != 1 && S/bc.sigma != 2 && S/bc.sigma != 4)
                throw std::runtime_error("S != alpha*sigma, where alpha=1, 2, or 4");
        }

    return 0;
}
