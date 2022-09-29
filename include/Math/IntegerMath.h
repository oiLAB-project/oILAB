/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_IntegerMath_h_
#define gbLAB_IntegerMath_h_

#include <numeric>
#include <Eigen/Dense>
#include <vector>
#include <iostream>


namespace gbLAB
{
    template <typename IntScalarType>
    struct IntegerMath
    {
        
        static IntScalarType sgn(const IntScalarType& a)
        {
            return a<0? -1 : 1;
        }
        
        static IntScalarType gcd(const IntScalarType &a, const IntScalarType &b)
        {
            const IntScalarType absA(abs(a));
            const IntScalarType absB(abs(b));
            return absB > 0 ? gcd(absB, absA % absB) : (absA > 0 ? absA : 1);
        }
        
        template<typename ArrayType>
        static IntScalarType gcd(const ArrayType& a)
        {
            switch (a.size())
            {
                case 0:
                {
                    throw std::runtime_error("gcd: array size is zero\n");
                    return 0;
                    break;
                }
                    
                case 1:
                {
                    return a[0];
                    break;
                }
                    
                case 2:
                {
                    return gcd(a[0],a[1]);
                    break;
                }
                    
                default:
                {
                    IntScalarType temp(a[0]);
                    for(long k=1;k<a.size();++k)
                    {
                        if (temp==0 && a[k]==0 && k!=a.size()-1)
                            temp= 0;
                        else
                            temp=gcd(temp,a[k]);
                    }
                    return temp;
                    break;
                }
            }
        }
        
        static IntScalarType lcm(const IntScalarType &a, const IntScalarType &b)
        {
            return a * b / gcd(a, b);
        }
        

        static Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic>
        integerGramSchmidt(const Eigen::Vector<IntScalarType,Eigen::Dynamic>& a)
        {
            int dim= a.size();
            Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> matrixA;
            Eigen::Matrix<IntScalarType,Eigen::Dynamic,1> tempx;
            matrixA.resize(1,1);
            matrixA << a(0);

            for(int i=0;i<dim-1;i++) {
                matrixA.conservativeResize(i + 1, i + 2);
                // zero the new column
                matrixA.col(i + 1) = Eigen::Vector<IntScalarType, Eigen::Dynamic>::Zero(i + 1);

                matrixA(0, i + 1) = a(i + 1);
                tempx.conservativeResize(i + 2, 1);

                // fill the new row with tempx
                if (i != 0)
                {
                    tempx(i + 1) = 0;
                    matrixA.row(i) = tempx;
                }

                Eigen::Vector<IntScalarType,Eigen::Dynamic> det=
                        Eigen::Vector<IntScalarType,Eigen::Dynamic>::Zero(i+2);
                for(int j=0;j<i+2; j++)
                {
                    Eigen::MatrixXd minor= Eigen::MatrixXd::Zero(i+1,i+1);
                    std::vector<IntScalarType> indicesToKeep(i+2);
                    std::iota(std::begin(indicesToKeep),std::end(indicesToKeep),0);
                    indicesToKeep.erase(std::remove(indicesToKeep.begin(),indicesToKeep.end(),j),indicesToKeep.end());

                    Eigen::Vector<IntScalarType,Eigen::Dynamic> indicesToKeepVector;
                    indicesToKeepVector= Eigen::Map<Eigen::Vector<IntScalarType,Eigen::Dynamic>>(indicesToKeep.data(),i+1);
                    minor= matrixA(Eigen::all, indicesToKeepVector).template cast<double>();
                    det(j)= round(minor.determinant());
                }

                IntScalarType e= gcd(det);


                if (det == Eigen::Vector<IntScalarType,Eigen::Dynamic>::Zero(i+2))
                {
                    tempx= Eigen::Vector<IntScalarType,Eigen::Dynamic>::Zero(i+2);
                    tempx(i)= 1;
                    continue;
                }

                for (int j = 0; j < i + 2; j++)
                    tempx(j) = pow(-1, j + 1) * det(j) / e;
            }
            matrixA.conservativeResize(dim,dim);
            matrixA.row(dim-1)= tempx;

            return matrixA.block(1,0,dim-1,dim);
        }
    };
}
#endif
