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
                                  const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                                  const BicrystalLatticeVectors& bicrystalConfig) :
          GbContinuum<dim>(getMesoStateGbDomain(mesoStateCslVectors),
                           get_xuPairs(gb,bs),
                           discretize(mesoStateCslVectors,gb),
                           bicrystalCoordsMap(gb,bicrystalConfig)),
          gb(gb),
          axis(axis),
          mesoStateCslVectors(mesoStateCslVectors),
          bicrystalConfig(bicrystalConfig),
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
    Eigen::Matrix<double, dim,dim-1> GbMesoState<dim>::getMesoStateGbDomain(const std::vector<LatticeVector<dim>>& mesoStateCslVectors)
    {
        Eigen::Matrix<double, dim,dim-1> mesoStateGbDomain;
        for(int i=1; i<dim; ++i)
            mesoStateGbDomain.col(i-1)= mesoStateCslVectors[i].cartesian();
        return mesoStateGbDomain;
    }

    template<int dim>
    std::map<OrderedTuplet<dim>, typename GbMesoState<dim>::VectorDimD>
          GbMesoState<dim>::get_xuPairs(const Gb<dim>& gb,
                                        const std::deque<std::pair<LatticeVector<dim>,VectorDimD>>& bs)
    {
        auto normal= gb.nA.cartesian();
        std::map<OrderedTuplet<dim>,VectorDimD> xuPairs;
        for(const auto& pair : bs)
        {
            OrderedTuplet<dim> key;
            VectorDimD value;
            // value = b/2
            value= pair.first.cartesian()/2;
            if(value.dot(normal)>0) {
                // replace value with a scaled normal inside this for loop
                // as a consequence, we will have to change OrderedTuple<dim> to OrderedTuple<Rational,dim> for this to generalize to ATGBs
                VectorDimD temp= pair.second - value;
                try {
                    //key << gb.bc.dscl.latticeVector(temp);
                    key << gb.bc.getLatticeVectorInD(gb.bc.A.latticeVector(temp));
                }
                catch(std::runtime_error& e)
                {
                    std::cout << e.what() << std::endl;
                    std::cout << "x = " << key.transpose() << ";   " << temp.transpose() << std::endl;
                    std::cout << "b = " << pair.first.cartesian().transpose()  << " ; s= " << pair.second.transpose() << std::endl;
                    exit(0);
                }
            }
            else {
                VectorDimD temp= pair.second + value;
                try {
                    //key << gb.bc.dscl.latticeVector(temp);
                    key << gb.bc.getLatticeVectorInD(gb.bc.B.latticeVector(temp));
                }
                catch(std::runtime_error& e)
                {
                    std::cout << e.what() << std::endl;
                    std::cout << "x = " << key.transpose() << ";   " << temp.transpose() << std::endl;
                    std::cout << "b = " << pair.first.cartesian().transpose()  << " ; s= " << pair.second.transpose() << std::endl;
                    exit(0);
                }
            }
            xuPairs[key]=value;
        }
        return xuPairs;
    }

    template<int dim>
    std::array<Eigen::Index,dim-1> GbMesoState<dim>::discretize(const std::vector<LatticeVector<dim>>& mesoStateCslVectors, const Gb<dim>& gb)
    {
        std::array<Eigen::Index,dim-1> n{};
        for(int i=1; i<dim; ++i)
            n[i-1]= 2*IntegerMath<int>::gcd(gb.bc.getLatticeVectorInD(mesoStateCslVectors[i]));
        return n;

    }


    template<int dim>
    std::map<OrderedTuplet<dim>,typename GbMesoState<dim>::VectorDimD> GbMesoState<dim>::bicrystalCoordsMap(const Gb<dim>& gb, const BicrystalLatticeVectors& bicrystalConfig)
    {
        std::map<OrderedTuplet<dim>,VectorDimD> idCoordsMap;

        for(const auto& latticeVector : bicrystalConfig) {
            auto latticeVectorInD= gb.bc.getLatticeVectorInD(latticeVector);
            OrderedTuplet<dim> key;
            key << latticeVectorInD;
            idCoordsMap[key]= latticeVectorInD.cartesian();
        }
        return idCoordsMap;

    }

    template<int dim>
    //template<int dm=dim>
    typename std::enable_if<dim==3,void>::type
    GbMesoState<dim>::box(const std::string& name) const
    {
        const auto& config= this->bicrystalConfig;
        std::vector<LatticeVector<3>> boxVectors;
        boxVectors.push_back(this->mesoStateCslVectors[0]);
        boxVectors.push_back(this->mesoStateCslVectors[1]);
        boxVectors.push_back(this->mesoStateCslVectors[2]);
        double boxHalfHeight= abs(boxVectors[0].cartesian().dot(gb.nA.cartesian().normalized()));

        std::vector<VectorDimD> referenceConfigA, deformedConfigA;
        std::vector<VectorDimD> referenceConfigB, deformedConfigB;
        std::vector<VectorDimD> configDscl;

        for (const auto &latticeVector: config) {
            if (&(latticeVector.lattice) == &(this->gb.bc.A) || &(latticeVector.lattice) == &(this->gb.bc.B))
            {
                VectorDimD x;
                OrderedTuplet<dim> temp;
                temp <<  gb.bc.getLatticeVectorInD(latticeVector);
                if (&(latticeVector.lattice) == &(this->gb.bc.A)) {
                    double height= latticeVector.cartesian().dot(gb.nA.cartesian().normalized());
                    if (height > 0 && height < boxHalfHeight-FLT_EPSILON)
                        temp << gb.bc.getLatticeVectorInD(gb.bc.B.latticeVector(latticeVector.cartesian() - 2 * height * gb.nA.cartesian().normalized()));
                }
                else if (&(latticeVector.lattice) == &(this->gb.bc.B)) {
                        double height = latticeVector.cartesian().dot(gb.nB.cartesian().normalized());
                        if (height > 0 && height < boxHalfHeight-FLT_EPSILON)
                            temp << gb.bc.getLatticeVectorInD(gb.bc.A.latticeVector(latticeVector.cartesian() - 2 * height * gb.nB.cartesian().normalized()));
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
