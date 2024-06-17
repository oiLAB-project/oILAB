//
// Created by Nikhil Chandra Admal on 5/27/24.
//

#include <GbMesoState.h>
#include <iostream>

namespace gbLAB {
    template<int dim>
    GbMesoState<dim>::GbMesoState(const Gb<dim>& gb,
                                  const ReciprocalLatticeVector<dim>& axis,
                                  const std::deque<std::pair<LatticeVector<dim>,VectorDimD>>& bs,
                                  const std::vector<LatticeVector<dim>>& gbCslVectors) :
                  GbContinuum<dim>(getGbDomain(gbCslVectors), get_xuPairs(bs,gb.nA.cartesian()), discretize(gbCslVectors,gb)),
                  gb(gb),
                  axis(axis),
                  gbCslVectors(gbCslVectors),
                  bs(bs)
    {
        // check
        /*
        auto xuPairs= get_xuPairs(bs);
        for(const auto& pair : xuPairs)
            std::cout << (pair.second-this->displacement(pair.first)).norm() << std::endl;
            */
    }

    template<int dim>
    Eigen::Matrix<double, dim,dim-1> GbMesoState<dim>::getGbDomain(const std::vector<LatticeVector<dim>>& gbCslVectors)
    {
        Eigen::Matrix<double, dim,dim-1> gbDomain;
        for(int i=0; i<dim-1; ++i)
            gbDomain.col(i)= gbCslVectors[i].cartesian();
        return gbDomain;
    }

    template<int dim>
    std::deque<std::pair<typename GbMesoState<dim>::VectorDimD, typename GbMesoState<dim>::VectorDimD>>
          GbMesoState<dim>::get_xuPairs(const std::deque<std::pair<LatticeVector<dim>,VectorDimD>>& bs,
                                        const VectorDimD& normal)
    {
              std::deque<std::pair<VectorDimD,VectorDimD>> xuPairs;
              std::pair<VectorDimD,VectorDimD> xu;
              for(const auto& pair : bs)
              {
                  xu.second= pair.first.cartesian()/2;
                  if(xu.second.dot(normal)>0)
                      xu.first= pair.second - xu.second;
                  else
                      xu.first= pair.second + xu.second;
                  xuPairs.push_back(xu);
              }
              return xuPairs;
    }

    template<int dim>
    std::array<Eigen::Index,dim-1> GbMesoState<dim>::discretize(const std::vector<LatticeVector<dim>>& gbCslVectors, const Gb<dim>& gb)
    {
        std::array<Eigen::Index,dim-1> n;
        for(int i=0; i<dim-1; ++i)
            n[i]= 2*IntegerMath<int>::gcd(gb.bc.getLatticeVectorInD(gbCslVectors[i]));
        //    n[1]= 2*IntegerMath<int>::gcd(gb.bc.getLatticeVectorInD(gbCslVectors[0]));

        return n;

    }

    template<int dim>
    //template<int dm=dim>
    typename std::enable_if<dim==3,void>::type
    GbMesoState<dim>::box(const int& heightFactor,
                        const int& dsclFactor,
                        const std::string& name) const
    {
        auto basis = this->gb.bc.csl.planeParallelLatticeBasis(
                this->gb.bc.getReciprocalLatticeDirectionInC(this->gb.nA.reciprocalLatticeVector()), true);
        LatticeVector<3> nonParallelC(basis[0].latticeVector());
        LatticeVector<3> vectorAlongAxisA(this->gb.bc.A.latticeDirection(this->axis.cartesian()).latticeVector());
        LatticeVector<3> vectorAlongAxisC(this->gb.bc.getLatticeDirectionInC(vectorAlongAxisA).latticeVector());

        std::vector<LatticeVector<3>> boxVectors;
        boxVectors.push_back(heightFactor * nonParallelC);
        boxVectors.push_back(gb.getPeriodVector(this->axis));
        boxVectors.push_back(vectorAlongAxisC);
        assert(abs(boxVectors[1].cartesian().dot(boxVectors[2].cartesian())) < FLT_EPSILON);

        auto config= this->gb.bc.box(boxVectors,1,1);
        std::vector<VectorDimD> referenceConfigA, deformedConfigA;
        std::vector<VectorDimD> referenceConfigB, deformedConfigB;
        std::vector<VectorDimD> configDscl;

        for (const auto &latticeVector: config) {
            if (&(latticeVector.lattice) == &(this->gb.bc.A) || &(latticeVector.lattice) == &(this->gb.bc.B))
            {
                VectorDimD x,temp;
                temp= latticeVector.cartesian();
                if (&(latticeVector.lattice) == &(this->gb.bc.A)) {
                    int planeIndex= latticeVector.dot(gb.nA);
                    double height;
                    if (planeIndex == 0)
                        height= FLT_EPSILON;
                    else
                        height= latticeVector.cartesian().dot(gb.nA.cartesian().normalized());
                    if (&(latticeVector.lattice) == &(this->gb.bc.A) && height > 0)
                        temp= latticeVector.cartesian() - 2 * height * gb.nA.cartesian().normalized();
                }
                else if (&(latticeVector.lattice) == &(this->gb.bc.B)) {
                        int planeIndex= latticeVector.dot(gb.nB);
                        double height;
                        if (planeIndex == 0)
                            height= FLT_EPSILON;
                        else
                            height = latticeVector.cartesian().dot(gb.nB.cartesian().normalized());
                        if (&(latticeVector.lattice) == &(this->gb.bc.B) && height > 0)
                            temp = latticeVector.cartesian() - 2 * height * gb.nB.cartesian().normalized();

                    }

                //x = latticeVector.cartesian() + this->displacement(latticeVector.cartesian());
                x = latticeVector.cartesian() + this->displacement(temp);

                if (&(latticeVector.lattice) == &(this->gb.bc.A) && x.dot(this->gb.nA.cartesian().normalized()) <= 1e-5)
                {
                    referenceConfigA.push_back(latticeVector.cartesian());
                    deformedConfigA.push_back(x);
                }
                else if (&(latticeVector.lattice) == &(this->gb.bc.B) && x.dot(this->gb.nB.cartesian().normalized()) <= 1e-5)
                    //else if (&(latticeVector.lattice) == &(gb.bc.B))
                {
                    referenceConfigB.push_back(latticeVector.cartesian());
                    deformedConfigB.push_back(x);
                }
            }

            // insert dscl atoms
            if (&(latticeVector.lattice) == &(this->gb.bc.dscl))
                configDscl.push_back(latticeVector.cartesian());
        }


        int nAtoms= referenceConfigA.size()+referenceConfigB.size()+configDscl.size();

        std::string referenceFile= name + "_reference0.txt";
        std::string deformedFile= name + "_reference1.txt";

        std::ofstream reference, deformed;
        reference.open(referenceFile);
        deformed.open(deformedFile);
        if (!reference || !deformed) std::cerr << "Unable to open files";
        reference << nAtoms << std::endl; deformed << nAtoms << std::endl;
        reference << "Lattice=\""; deformed << "Lattice=\"";

        reference << std::setprecision(15) << (2*boxVectors[0].cartesian()).transpose() << " ";
        deformed << std::setprecision(15) << (2*boxVectors[0].cartesian()).transpose() << " ";
        reference << std::setprecision(15) << (boxVectors[1].cartesian()).transpose() << " ";
        deformed << std::setprecision(15) << (boxVectors[1].cartesian()).transpose() << " ";
        reference << std::setprecision(15) << (boxVectors[2].cartesian()).transpose();
        deformed << std::setprecision(15) << (boxVectors[2].cartesian()).transpose();
        reference << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\""; deformed << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
        reference << std::setprecision(15) << (-1 * boxVectors[0].cartesian()).transpose() << "\"" << std::endl;
        deformed << std::setprecision(15) << (-1 * boxVectors[0].cartesian()).transpose() << "\"" << std::endl;

        for(const auto& position : referenceConfigA)
            reference << 1 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : referenceConfigB)
            reference << 2 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : deformedConfigA)
            deformed << 1 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : deformedConfigB)
            deformed << 2 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : configDscl) {
            reference << 4 << " " << std::setprecision(15) << position.transpose() << "  " << 0.01 << std::endl;
            deformed << 4 << " " << std::setprecision(15) << position.transpose() << "  " << 0.01 << std::endl;
        }

        reference.close();
        deformed.close();
    }

    //template class GbMesoState<2>;
    template class GbMesoState<3>;
}
