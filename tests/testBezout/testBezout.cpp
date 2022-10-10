#include <IntegerMath.h>
#include <random>


using namespace gbLAB;
int main()
{
    typedef long long int IntScalarType;

    int dim=4;
    for (int trial=0; trial<100; trial++) {
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

        std::cout << "Output = " << out.transpose() << std::endl;
        std::cout << "dot product = " << out.dot(in) << std::endl;

        try {
            if (out.dot(in) != 1)
                throw std::runtime_error("Test Bezout failed");
        }
        catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return -1;
        }
    }
    return 0;

}
