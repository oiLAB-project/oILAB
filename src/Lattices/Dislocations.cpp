//
// Created by Nikhil Chandra Admal on 2/18/24.
//
#include <Dislocations.h>
#include<cfloat>
#include <iostream>
#include <randomInteger.h>

namespace gbLAB {

    Dislocations::Dislocations(const double& a2,
                               const double& shiftSize,
                               const int& nImages):
    a2(a2),
    shift((Vector2d()<< shiftSize,0).finished()),
    nImages(nImages){}

    void Dislocations::insertDislocationDipole(Matrix2d &ends, Vector2d &b)
    {
        (*this).emplace_back(DislocationDipole(ends,b,shift,nImages));
    }

    void Dislocations::removeDislocationDipole(const int& index)
    {
        this->erase(this->begin() + index);
    }
    const std::vector<DislocationDipole>& Dislocations::dislocations() const
    {
        return *this;
    }

    // Perform a line search along the Newton direction
    double Dislocations::line_search(const Vector2d &x, const Vector2d &dx, const Vector2d &X, const int& branch) const {
        double alpha = 1.0;
        double beta = 0.5; // Backtracking factor
        double c = 0.01; // Sufficient decrease parameter

        double fx_norm = (inverseDeformationMap(x,branch) - X).norm();
        Vector2d x_new = x + alpha * dx;
        double fx_new_norm = (inverseDeformationMap(x_new,branch) - X).norm();
        //while (fx_new_norm > (1 - c * alpha) * fx_norm && alpha > FLT_EPSILON) {
        int iter= 0;
        const int max_iterations = 6; // Maximum number of iterations
        while ((fx_new_norm >= (1 - c * alpha) * fx_norm &&
               iter<max_iterations) || x_new(1)*X(1)<0  )
        {
            alpha *= beta;
            x_new = x + alpha * dx;
            fx_new_norm = (inverseDeformationMap(x_new, branch) - X).norm();
            iter++;
        }

        return alpha;
    }

    typename Dislocations::Vector2d Dislocations::deformationMap(const Vector2d &X, const int& branch) const {
        const double tolerance = 1e-5;
        const int max_iterations = 1000;

        Vector2d x = X;
        int iterations = 0;
        double error = tolerance + 1;

        for(const auto& dipole : dislocations()) {
            for (int i = -dipole.nImages; i <= dipole.nImages; ++i) {
                const auto& A(Vector2d(dipole.ends.col(0)));
                const auto& B(Vector2d(dipole.ends.col(1)));

                const auto lp(dipole.localPosition(A, B, x + i * dipole.shift));
                const Vector2d &xL(lp.first);
                if (abs(xL(1)) < FLT_EPSILON) {
                    return x;
                }
            }
        }

        while (error > tolerance && iterations < max_iterations) {
            Vector2d fx = inverseDeformationMap(x,branch);
            Matrix2d J = jacobian(x,branch);

            // Solve J*dx = y - f(x) for dx using LU decomposition
            Vector2d dx = J.lu().solve(X - fx);

            // Perform line search along the Newton direction
            double alpha = line_search(x, dx, X, branch);


            // Update x
            x += alpha * dx;

            // Calculate the error
            //error = dx.norm();
            error= (X-inverseDeformationMap(x,branch)).norm();
            iterations++;
        }

        if (iterations == max_iterations) {
            /*
            std::cout << "branch = " << branch << std::endl;
            std::cout << "X = " << X.transpose() << std::endl;
            std::cout << "x = " << x.transpose() << std::endl;
            std::cout << "fx = " << inverseDeformationMap(x,branch).transpose() << std::endl;
            std::cout << "Error = " << error << std::endl;
             */
            throw std::runtime_error("Newton's method did not converge within the maximum number of iterations.");
        }

        return x;
    }

    typename Dislocations::Matrix2d Dislocations::jacobian(const Vector2d &x, const int& branch) const {
        Matrix2d temp(Matrix2d::Zero());

        for(const auto& dipole : dislocations())
            temp+= dipole.b * dipole.solidAngleGradient(x,branch).transpose() / (4 * M_PI);
        return Eigen::Matrix2d::Identity() + temp;
    }

    typename Dislocations::Vector2d Dislocations::inverseDeformationMap(const Vector2d& x, const int& branch) const
    {
        Vector2d temp(Vector2d::Zero());
        for(const auto& dipole : dislocations())
            temp= temp + dipole.displacement(x,branch);
        return x - temp;
    }

    double Dislocations::edgeElasticEnergyKernel(const Vector2d& x1,const Vector2d& x2,const Vector2d& b1,const Vector2d& b2) const
    {
        const Vector2d r(x2-x1);
        const double ra2(r.squaredNorm()+a2);
        const double ra(sqrt(ra2));
        const double a(sqrt(a2));
        return -(log(ra/a)-a2/ra2)*b1.dot(b2)+b1.dot(r)*b2.dot(r)/ra2;
    }
    double Dislocations::elasticEnergy() const
    {
        double temp(0.0);


        // Collect end points and Burgers vectors
        std::vector<std::pair<Vector2d,Vector2d>> singleDislocations; // position and Burgers
        for(const auto& dipole : dislocations())
        {
            singleDislocations.emplace_back(dipole.ends.col(0), dipole.b); // "positive"
            singleDislocations.emplace_back(dipole.ends.col(1),-dipole.b); // "negative"
        }

        for(const auto& dis1 : singleDislocations)
        {
            for(const auto& dis2 : singleDislocations)
            {
                for(int i=-nImages;i<=nImages;++i)
                {
                    temp+=edgeElasticEnergyKernel(dis1.first,dis2.first+i*shift,dis1.second,dis2.second);
                }
            }
        }

        //return mu/(4.0*M_PI*(1-nu))*temp;
        return temp;
    }





}