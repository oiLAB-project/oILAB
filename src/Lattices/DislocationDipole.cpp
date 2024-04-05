//
// Created by Nikhil Chandra Admal on 2/17/24.
//
#include<DislocationDipole.h>
#include<cfloat>
namespace gbLAB
{
    DislocationDipole::DislocationDipole(const Matrix2d& ends_in,
                                         const Vector2d& b_in,
                                         const Vector2d& shift_in,
                                         const int& nImages_in) :
                                         ends(ends_in),
                                         b(b_in),
                                         shift(shift_in),
                                         nImages(nImages_in){}

    std::pair<typename DislocationDipole::Vector2d,double> DislocationDipole::localPosition(const Vector2d& A,const Vector2d& B,const Vector2d& x)
    {
        const Vector2d A2B(B-A);
        const double normA2B(A2B.norm());
        const Vector2d t(A2B/normA2B);
        const Vector2d n((Vector2d()<<-t(1),t(0)).finished());
        const Vector2d c(0.5*(A+B));
        return std::make_pair((Vector2d()<<(x-c).dot(t),(x-c).dot(n)).finished(),0.5*normA2B);
    }

    double DislocationDipole::solidAngle(const Vector2d& x, const int& branch) const
    {
        double temp(0.0);
        for(int i=-nImages;i<=nImages;++i)
        {
            const auto& A(Vector2d(ends.col(0)));
            const auto& B(Vector2d(ends.col(1)));

            const auto lp(localPosition(A,B,x+i*shift));
            const Vector2d& xL(lp.first);
            const double& halfLength(lp.second);

            if(abs(xL(1))<FLT_EPSILON || branch < 0)
            //if(branch < 0)
            {
                if(abs(xL(0))<halfLength)
                {
                    if (xL(1)>0 && branch > 0|| (xL(1)<0 && branch < 0))
                        temp+=2.0*M_PI;
                    else if (xL(1)<0 && branch > 0|| (xL(1)>0 && branch < 0))
                        temp+=-2.0*M_PI;
                    else
                        temp+=0.0;
                    //temp+=2.0*M_PI;
                }
            }
            else
            {
                const double Yterm(abs(1.0/xL(1)));
                //temp+=2.0*sgn(xL(1))*(std::atan((xL(0)-halfLength)*Yterm)-std::atan((xL(0)+halfLength)*Yterm));
                temp+=2.0*sgn(xL(1))*(-std::atan((xL(0)-halfLength)*Yterm)+std::atan((xL(0)+halfLength)*Yterm));

            }
        }
        return temp;
    }

    typename DislocationDipole::Vector2d DislocationDipole::displacement(const Vector2d & x, const int& branch) const
    {
        return -(solidAngle(x,branch)/4.0/M_PI)*b;
    }

    typename DislocationDipole::Vector2d DislocationDipole::solidAngleGradient(const Vector2d& x, const int& branch) const
    {
        Vector2d temp(Vector2d::Zero());
        for(int i=-nImages;i<=nImages;++i)
        {
            const auto& A(Vector2d(ends.col(0)));
            const auto& B(Vector2d(ends.col(1)));

            const auto lp(localPosition(A,B,x+i*shift));
            const Vector2d& xL(lp.first);
            const double& halfLength(lp.second);

            if(abs(xL(1))<FLT_EPSILON || branch<0)
            //if(branch<0)
            {
                if(abs(xL(0))<halfLength)
                {
                    temp+= Vector2d(0,0);
                }
            }
            else
            {
                const double Yterm(1.0/xL(1));
                const double Xterm1(xL(0)-halfLength);
                const double Xterm2(xL(0)+halfLength);
                //const double factor= 2.0* (-1.0/(1+std::pow(Xterm1*Yterm,2)) +1.0/(1+std::pow(Xterm2*Yterm,2)));
                const double factor= 2.0* (
                       -1.0/(std::pow(Xterm1,2) + std::pow(xL(1),2)) +
                        1.0/(std::pow(Xterm2,2) + std::pow(xL(1),2))
                                          );
                const Vector2d A2B(B-A);
                const double normA2B(A2B.norm());
                const Vector2d t(A2B/normA2B);
                const Vector2d n((Vector2d()<<-t(1),t(0)).finished());
                //temp+= factor*(t*Yterm -xL(0)*std::pow(Yterm,2)*n);
                temp+= factor*(t*xL(1) -xL(0)*n);
                //temp+=2.0*sgn(xL(1))*(-std::atan((xL(0)-halfLength)*Yterm)+std::atan((xL(0)+halfLength)*Yterm));
            }

        }
        return temp;

    }

}

