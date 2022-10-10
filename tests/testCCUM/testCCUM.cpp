#include <IntegerMath.h>
#include <random>


using namespace gbLAB;
int main()
{
    typedef long long int IntScalarType;

    int dim=14;
    for (int trial=0; trial<1000; trial++) {
        Eigen::Vector<IntScalarType, Eigen::Dynamic> in(dim);

        // generate random input direction
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); // seed the generator
        std::uniform_int_distribution<> distr(-100, 100); // define the range
        for (int n = 0; n < dim; ++n) {
            in(n) = (IntScalarType) distr(gen);
        }
        in = in / IntegerMath<IntScalarType>::gcd(in);

        std::cout << "Input = " << in.transpose() << std::endl;
        auto out = IntegerMath<IntScalarType>::solveBezout(in);
        std::cout << "Bezout solver output = " << out.transpose() << std::endl;

        auto matrix = IntegerMath<IntScalarType>::ccum(out);
        int determinant = round(matrix.template cast<double>().determinant());
        int dot_product= matrix.col(0).dot(in);

        std::cout << "CCUM solver output = " << std::endl;
        std::cout << matrix << std::endl;
        std::cout << "dot product = " << dot_product << std::endl;
        std::cout << "determinant = " << determinant << std::endl;
        try {
            if (abs(dot_product) != 1 || abs(determinant) != 1)
            {
                throw std::runtime_error("Test CCUM failed");
            }
        }

        catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return -1;
        }
    }
    return 0;

}
