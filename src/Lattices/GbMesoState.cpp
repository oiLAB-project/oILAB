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
                           get_xuPairs(gb,mesoStateCslVectors,bs),
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
    std::map<OrderedTuplet<dim+1>, typename GbMesoState<dim>::VectorDimD>
          GbMesoState<dim>::get_xuPairs(const Gb<dim>& gb,
                                        const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                                        const std::deque<std::pair<LatticeVector<dim>,VectorDimD>>& bs)
    {
        auto normal= gb.nA.cartesian();
        std::map<OrderedTuplet<dim+1>,VectorDimD> xuPairs;
        std::vector<LatticeVector<dim>> bicrystalBoxVectors(mesoStateCslVectors);
        bicrystalBoxVectors[0]= 2*mesoStateCslVectors[0];
        VectorDimD shift;
        shift << -0.5-FLT_EPSILON,-FLT_EPSILON,-FLT_EPSILON;


        for(const auto& [b,s] : bs)
        {
            OrderedTuplet<dim+1> keyx;
            VectorDimD valueu;
            // value = b/2
            valueu= b.cartesian()/2;
            VectorDimD tempx= s-valueu;


            if(tempx.dot(normal)<0)
            {
                // modulo tempx w.r.t the bicrystal box
                LatticeVector<dim>::modulo(tempx,bicrystalBoxVectors,shift);
                try {
                    keyx << gb.bc.getLatticeVectorInD(gb.bc.A.latticeVector(tempx)),1;
                }
                catch(std::runtime_error& e)
                {
                    std::cout << e.what() << std::endl;
                    std::cout << "x key = " << keyx.transpose() << ";   " << tempx.transpose() << std::endl;
                    std::cout << "b = " << b.cartesian().transpose()  << " ; s= " << s.transpose() << std::endl;
                    exit(0);
                }

            }
            else
            {
                tempx= tempx + 2*valueu.dot(gb.nA.cartesian().normalized()) * gb.nA.cartesian().normalized();
                // modulo tempx w.r.t the bicrystal box
                LatticeVector<dim>::modulo(tempx,bicrystalBoxVectors,shift);

                try {
                    keyx << gb.bc.getLatticeVectorInD(gb.bc.B.latticeVector(tempx)),2;
                }
                catch(std::runtime_error& e)
                {
                    std::cout << e.what() << std::endl;
                    std::cout << "x key = " << keyx.transpose() << ";   " << tempx.transpose() << std::endl;
                    std::cout << "b = " << b.cartesian().transpose()  << " ; s= " << s.transpose() << std::endl;
                    exit(0);
                }

            }
            xuPairs[keyx]=valueu;
        }
        return xuPairs;
    }

    template<int dim>
    std::array<Eigen::Index,dim-1> GbMesoState<dim>::discretize(const std::vector<LatticeVector<dim>>& mesoStateCslVectors, const Gb<dim>& gb)
    {
        std::array<Eigen::Index,dim-1> n{};
        for(int i=1; i<dim; ++i)
            n[i-1]= 4*IntegerMath<int>::gcd(gb.bc.getLatticeVectorInD(mesoStateCslVectors[i]));
        return n;

    }


    template<int dim>
    std::map<OrderedTuplet<dim+1>,typename GbMesoState<dim>::VectorDimD> GbMesoState<dim>::bicrystalCoordsMap(const Gb<dim>& gb, const BicrystalLatticeVectors& bicrystalConfig)
    {
        std::map<OrderedTuplet<dim+1>,VectorDimD> idCoordsMap;

        for(const auto& latticeVector : bicrystalConfig) {
            auto latticeVectorInD= gb.bc.getLatticeVectorInD(latticeVector);
            OrderedTuplet<dim+1> key;
            if (&latticeVector.lattice == &gb.bc.A)
                key << latticeVectorInD,1;
            else if (&latticeVector.lattice == &gb.bc.B)
                key << latticeVectorInD,2;
            idCoordsMap[key]= latticeVectorInD.cartesian();
        }
        return idCoordsMap;

    }

    /*-------------------------------------*/
    template<int dim>
    double GbMesoState<dim>::energy() const
    {
        // form box
        // run a python script to calculate energy
        return 0.0;
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

        std::vector<VectorDimD> referenceConfigA, deformedConfigA;
        std::vector<VectorDimD> referenceConfigB, deformedConfigB;
        std::vector<VectorDimD> configDscl;

        double bmax= 0.0;
        for(const auto& elem : bs)
            bmax= max(bmax,elem.first.cartesian().norm());

        for (const auto &latticeVector: config) {
            VectorDimD x;
            OrderedTuplet<dim+1> temp;
            if (&(latticeVector.lattice) == &(this->gb.bc.A)) {
                temp <<  gb.bc.getLatticeVectorInD(latticeVector),1;
                double height= latticeVector.cartesian().dot(gb.nA.cartesian().normalized());
                LatticeDirection<dim> latticeVectorInAalongnA(gb.bc.A.latticeDirection(gb.nA.cartesian()));
                LatticeDirection<dim> dsclDirectionAlongnA(gb.bc.getLatticeDirectionInD(latticeVectorInAalongnA.latticeVector()));
                LatticeVector<dim> dsclVector= round(abs(2.0*height/dsclDirectionAlongnA.cartesian().norm())) * dsclDirectionAlongnA.latticeVector();

                if (height > 0 && height < bmax) // latticeVector is outside A
                    //temp << gb.bc.getLatticeVectorInD(gb.bc.B.latticeVector(latticeVector.cartesian() - 2.0 * height * gb.nA.cartesian().normalized())),2;
                    temp << gb.bc.getLatticeVectorInD(latticeVector)-dsclVector, 2;
            }
            else if (&(latticeVector.lattice) == &(this->gb.bc.B)) {
                temp <<  gb.bc.getLatticeVectorInD(latticeVector),2;
                double height = latticeVector.cartesian().dot(gb.nB.cartesian().normalized());
                LatticeDirection<dim> latticeVectorInBalongnB(gb.bc.B.latticeDirection(gb.nB.cartesian()));
                LatticeDirection<dim> dsclDirectionAlongnB(gb.bc.getLatticeDirectionInD(latticeVectorInBalongnB.latticeVector()));
                LatticeVector<dim> dsclVector= round(abs(2.0*height/dsclDirectionAlongnB.cartesian().norm())) * dsclDirectionAlongnB.latticeVector();
                if (height > 0 && height < bmax) // latticeVector is outside B
                    //temp << gb.bc.getLatticeVectorInD(gb.bc.A.latticeVector(latticeVector.cartesian() - 2.0 * height * gb.nB.cartesian().normalized())),1;
                    temp << gb.bc.getLatticeVectorInD(latticeVector)-dsclVector, 1;
            }

            //x = latticeVector.cartesian() + this->displacement(latticeVector.cartesian());

            if (&(latticeVector.lattice) == &(this->gb.bc.A)) {
                x = latticeVector.cartesian() + this->displacement(temp) + this->uAverage;
            }
            else if (&(latticeVector.lattice) == &(this->gb.bc.B) )
                x = latticeVector.cartesian() + this->displacement(temp) - this->uAverage;
            else
                x = latticeVector.cartesian() + this->displacement(temp);


            if (&(latticeVector.lattice) == &(this->gb.bc.A) && x.dot(this->gb.nA.cartesian().normalized()) <= FLT_EPSILON)
            //if (&(latticeVector.lattice) == &(this->gb.bc.A))
            {
                referenceConfigA.push_back(latticeVector.cartesian());
                deformedConfigA.push_back(x);
            }
            else if (&(latticeVector.lattice) == &(this->gb.bc.B) && x.dot(this->gb.nB.cartesian().normalized()) <= FLT_EPSILON)
            //else if (&(latticeVector.lattice) == &(this->gb.bc.B))
            {
                referenceConfigB.push_back(latticeVector.cartesian());
                deformedConfigB.push_back(x);
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
