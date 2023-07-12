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
        
        template<typename T>
        static IntScalarType gcd(const Eigen::MatrixBase<T>& a)
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
                    return a(0);
                    break;
                }
                    
                case 2:
                {
                    return gcd(a(0),a(1));
                    break;
                }
                    
                default:
                {
                    IntScalarType temp(a(0));
                    for(long k=1;k<a.size();++k)
                    {
                        if (temp==0 && a(k)==0 && k!=a.size()-1)
                            temp= 0;
                        else
                            temp=gcd(temp,a(k));
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

        template<typename T>
        static IntScalarType lcm(const Eigen::MatrixBase<T>& a)
        {
            switch (a.size())
            {
                case 0:
                {
                    throw std::runtime_error("lcm: array size is zero\n");
                    return 0;
                    break;
                }

                case 1:
                {
                    return a(0);
                    break;
                }

                case 2:
                {
                    return lcm(a(0),a(1));
                    break;
                }

                default:
                {
                    IntScalarType temp(a(0));
                    for(long k=1;k<a.size();++k)
                    {
                        if (temp==0 && a(k)==0 && k!=a.size()-1)
                            temp= 0;
                        else
                            temp=lcm(temp,a(k));
                    }
                    return temp;
                    break;
                }
            }

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
                    minor= matrixA(Eigen::indexing::all, indicesToKeepVector).template cast<double>();
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

        // Finds such x and y, that a  x + b  y = gcd(a, b)
        //https://github.com/ADJA/algos/blob/master/NumberTheory/DiophantineEquation.cpp
        static IntScalarType extended_gcd(IntScalarType a, IntScalarType b, IntScalarType &x, IntScalarType &y)
        {
            if (b == 0)
            {
                x = 1;
                y = 0;
                return a;
            }
            IntScalarType x1, y1;
            IntScalarType g = extended_gcd(b, a % b, x1, y1);
            x = y1;
            y = x1 - (a / b) * y1;
            return g;
        }

        // Solves equation a  x + b  y = c, writes answer to x and y
        static void solveDiophantine2vars(IntScalarType a, IntScalarType b, IntScalarType c, IntScalarType &x, IntScalarType &y)
        {

            IntScalarType g = extended_gcd(a, b, x, y);

            if (c % g != 0)
            {
                std::cout << a << "  " << b << "  " << c << "   " << g << std::endl;
                puts("Impossible");
                exit(0);
            }

            c /= g;

            x = x * c ;
            y = y * c ;
        }

        // Find u such that a_1 u_1 + a_2 u_2 + .... + a_n u_n = 1
        // Method 0:
        // "A fast algorithm to find reduced hyperplane unit cells and solve N-dimensional Bézout's identities."
        //  Acta Crystallographica Section A: Foundations and Advances 77.5 (2021).
        // template<typename ArrayType>
        // static ArrayType solveBezout(const ArrayType& a)
        template<typename T>
        static Eigen::Vector<IntScalarType,Eigen::Dynamic> solveBezout(const Eigen::MatrixBase<T>& a)
        {
            int n= a.size();
            if (n<2) throw std::runtime_error("the size of arrays should be at least two");
            if (abs(IntegerMath<IntScalarType>::gcd(a)) != 1 || a.isZero())
                throw std::runtime_error("No solution since the gcd is not 1 or -1.");

            Eigen::Vector<IntScalarType,Eigen::Dynamic> u(n);

            IntScalarType p,q;
            IntScalarType g12= extended_gcd(a(0),a(1),p,q); // g12 - gcd of a(0) and a(1)

            // form the n-1 sized vector na= {g,a2,...,a_{n-1}}
            Eigen::Vector<IntScalarType,Eigen::Dynamic> na(n-1);
            na(0)= g12; na(Eigen::seq(1,n-2))= a(Eigen::seq(2,n-1));

            Eigen::Vector<IntScalarType,Eigen::Dynamic> k(n-1);
            if (n==2)
            {
                IntScalarType g= extended_gcd(a(0),a(1),u(0),u(1));
                if (g<0) u= -u;
                return u;
            }
            else
            {
                k = solveBezout(na);
                // {p*k0,q*k0,k1,..,k_{n-2}}
                u(0) = p * k(0);
                u(1) = q * k(0);
                u(Eigen::seq(2, n - 1)) = k(Eigen::seq(1, n - 2));
                return u;
            }
        }

        template<typename T>
        static Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> ccum(const Eigen::MatrixBase<T>& qin)
        //static Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic>
        //ccum(const Eigen::Vector<IntScalarType,Eigen::Dynamic>& qin)
        {
            int n= qin.size();
            Eigen::Vector<IntScalarType,Eigen::Dynamic> q(n);
            q= qin;

            try
            {
                if (n==1)
                    throw std::runtime_error("Size of the input matrix for CCUM should be >1.");
                if (abs(gcd(qin)) != 1)
                    throw std::runtime_error("The absolute gcd of the input array is not 1.");
            }
            catch(std::runtime_error& e)
            {
                std::cerr << e.what();
            }

            Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> Q=
                    Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic>::Identity(n,n);

            // first ensure that q(0) is non-zero; otherwise swap
            bool swapFirst= false;
            int swapFirstElementWith;
            if (q(0)==0)
            {
                int elementIndex= -1;
                for (auto element : q)
                {
                    elementIndex++;
                    if (element !=0) {
                        // swap
                        IntScalarType temp= q(0);
                        q(0)= element;
                        q(elementIndex)= temp;
                        swapFirst = true;
                        swapFirstElementWith= elementIndex;
                        break;
                    }
                }
            }
            Q.col(0)= q;

            int elementIndex= -1;
            IntScalarType y;
            // ensure q(0) and q(1) are co-prime; otherwise swap
            for (const IntScalarType& element : q)
            {
                elementIndex++;
                if (elementIndex == 0) continue;
                if (!(q(0)==0 && element==0))
                    y= gcd(q(0),element);
                else
                    continue;
                if (abs(y)==1)
                {
                    // swap rows
                    int temp= Q(1,0);
                    Q(1,0)= Q(elementIndex,0);
                    Q(elementIndex,0)= temp;

                    IntScalarType u,v;
                    IntegerMath<IntScalarType>::solveDiophantine2vars(q(0),element,1,u,v);
                    //Q(0,elementIndex)=-v; Q(elementIndex,elementIndex)= u;
                    Q(0,1)=-v; Q(1,1)= u;
                    // unswap
                    Q.row(1).swap(Q.row(elementIndex));

                    if(swapFirst) Q.row(0).swap(Q.row(swapFirstElementWith));
                    return Q;
                }
            }

            // At this point, q(0) is not co-prime with any other elements of q
            if (!(q(0) ==0 && q(1) == 0))
                y= gcd(q(0),q(1));
            else
                y= 0;
            IntScalarType u,v;
            IntegerMath<IntScalarType>::solveDiophantine2vars(q(0),q(1),y,u,v);

            IntScalarType alpha,beta;
            IntegerMath<IntScalarType>::solveDiophantine2vars(u,v,1,alpha,beta);
            Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> M=
                    Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic>::Identity(n,n);
            M.block(0,0,2,2) << u,     v,
                               -beta, alpha;

            Eigen::Vector<IntScalarType,Eigen::Dynamic> t;
            t= M*q;
            elementIndex= -1;
            int swappedElementIndex;
            bool swapped= false;
            for (IntScalarType& element : t)
            {
                elementIndex++;
                if (elementIndex<2) continue;
                // Find a element in {t(2),...,} that y:=gcd(q(0), q(1)) does not divide
                // We are guaranteed to find such an element since gcd(q)=1
                if (element%y != 0)
                {
                    // swap the found element with t(1)
                    int temp= t(elementIndex);
                    t(elementIndex)= t(1);
                    t(1)= temp;
                    swappedElementIndex= elementIndex;
                    swapped= true;
                    break;
                }
                assert(elementIndex!=t.size());
            }

            Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> tempMatrix(n,n) ;
            tempMatrix= ccum(t);

            // inv(P)*ccum(t)
            // swap rows 0 and swappedElementIndex
            if (swapped) tempMatrix.row(1).swap(tempMatrix.row(swappedElementIndex));

            // inv(M)*inv(P)*ccum(t)
            Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> invM=
                    Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic>::Identity(n,n);
            invM.block(0,0,2,2) << alpha, -v,
                                   beta,   u;
            Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> output(n,n);
            output= invM*tempMatrix;
            if(swapFirst) output.row(0).swap(output.row(swapFirstElementWith));
            return output;
        }

    };

    template <typename IntScalarType, int dim>
    class MatrixDimIExt : public Eigen::Matrix<IntScalarType,dim,dim>
    {
    public:
        static Eigen::Matrix<IntScalarType,dim,dim> adjoint(const Eigen::Matrix<IntScalarType,dim,dim>& input)
        {
            Eigen::Matrix<IntScalarType,dim,dim> output;
            for (int i=0; i<dim; i++)
            {
                // remove i-th row
                Eigen::Matrix<IntScalarType,dim-1,dim> temp1 ;
                temp1 << input(Eigen::seq(0,i-1),Eigen::all),
                         input(Eigen::seq(i+1,dim-1),Eigen::all);
                for (int j=0; j<dim; j++)
                {
                    // remove jth column
                    Eigen::Matrix<IntScalarType,dim-1,dim-1> temp2 ;
                    temp2 << temp1(Eigen::all,Eigen::seq(0,j-1)), temp1(Eigen::all,Eigen::seq(j+1,dim-1));

                    output(i,j)=(IntScalarType) std::pow(-1,i+j+2)*temp2.template cast<double>().determinant();
                }
            }

            return output.transpose();
        }
    };
}
#endif
