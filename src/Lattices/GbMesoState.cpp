//
// Created by Nikhil Chandra Admal on 5/27/24.
//

#include<Lammps.h>
#include <GbMesoState.h>
#include <iostream>
#include <Python.h>

namespace gbLAB {
    template<int dim>
    GbMesoState<dim>::GbMesoState(const Gb<dim>& gb,
                                  const ReciprocalLatticeVector<dim>& axis,
                                  const std::deque<std::tuple<LatticeVector<dim>,VectorDimD,int>>& bs,
                                  const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                                  const BicrystalLatticeVectors& bicrystalConfig)
                                  try :
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
    catch(std::runtime_error& e)
    {
        //throw(std::runtime_error("GB Mesostate construction failed"));
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
                                        const std::deque<std::tuple<LatticeVector<dim>,VectorDimD,int>>& bs)
    {
        auto normal= gb.nA.cartesian().normalized();
        std::map<OrderedTuplet<dim+1>,VectorDimD> xuPairs;
        std::vector<LatticeVector<dim>> bicrystalBoxVectors(mesoStateCslVectors);
        bicrystalBoxVectors[0]= 2*mesoStateCslVectors[0];
        VectorDimD shift;
        shift << -0.5-FLT_EPSILON,-FLT_EPSILON,-FLT_EPSILON;


        for(const auto& [b,s,include] : bs)
        {
            OrderedTuplet<dim+1> keyx;
            VectorDimD valueu;
            // value = b/2
            valueu= b.cartesian()/2;
            VectorDimD tempx= s-valueu;


            /*
            if(tempx.dot(normal)<FLT_EPSILON) // tempx is in lattice A
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
            else    // tempx is in lattice B
            {
                tempx= tempx + 2*valueu.dot(normal) * normal;
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
             */
            // modulo tempx w.r.t the bicrystal box
            LatticeVector<dim>::modulo(tempx,bicrystalBoxVectors,shift);
            try {
                if (tempx.dot(normal) < FLT_EPSILON) // tempx is in lattice A
                    keyx << gb.bc.getLatticeVectorInD(gb.bc.A.latticeVector(tempx)), 1;
                else    // tempx is in lattice B region
                    keyx << gb.bc.getLatticeVectorInD(gb.bc.A.latticeVector(tempx)), -1;
            }
            catch(std::runtime_error& e) {
                std::cout << e.what() << std::endl;
                std::cout << "x key = " << keyx.transpose() << ";   " << tempx.transpose() << std::endl;
                std::cout << "b = " << b.cartesian().transpose()  << " ; s= " << s.transpose() << std::endl;
                exit(0);
            }
            xuPairs[keyx]=valueu;
        }
        if(xuPairs.size() != bs.size())
            throw(std::runtime_error("Clash in constraints."));
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
            if (&latticeVector.lattice == &gb.bc.A) {
                if (latticeVector.dot(gb.nA) <= 0)
                    key << latticeVectorInD, 1;
                else
                    key << latticeVectorInD, -1;
            }
            else if (&latticeVector.lattice == &gb.bc.B) {
                if (latticeVector.dot(gb.nB) <= 0)
                    key << latticeVectorInD, 2;
                else
                    key << latticeVectorInD, -2;
            }
            idCoordsMap[key]= latticeVectorInD.cartesian();
        }
        return idCoordsMap;

    }

    /*-------------------------------------*/
    template<int dim>
    std::pair<double,double> GbMesoState<dim>::densityEnergy() const
    {
        //std::string lammpsLocation= "/Users/Nikhil/Documents/Academic/Software/lammps-29Aug2024/build/lmp";
        std::string lammpsLocation= "/Users/Nikhil/Documents/Academic/Software/lammps-15May15/src/lmp_serial";
        //std::string lammpsLocation= "/usr/bin/lmp";
        box("temp" + std::to_string(omp_get_thread_num()));
        //auto xx= densityEnergyPython();
        auto yy= energy( lammpsLocation, "temp" + std::to_string(omp_get_thread_num()) + "_reference1.txt","Cu_mishin1.eam.alloy");
        //return energy( lammpsLocation, "temp_reference1.txt","Cu_mishin1.eam.alloy");
	return yy;
    }



    /*-------------------------------------*/
 template<int dim>
 std::pair<double,double> GbMesoState<dim>::densityEnergyPython() const
 {
     // form box
     // run a python script to calculate energy

     setenv("PYTHONPATH", ".", 1);
     if(!Py_IsInitialized()) {
         std::cout << "Initializing Python Interpreter" << std::endl;
         Py_Initialize();
     }

     box("temp");
     PyObject* pyModuleString = PyUnicode_FromString((char*)"lammps");
     PyObject* pyModule = PyImport_Import(pyModuleString);
     if (pyModule == nullptr)
     {
         PyErr_Print();
         std::cout << "cannot import python script" << std::endl;
         std::exit(0);
     }
     PyObject* pDict = PyModule_GetDict(pyModule);
     PyObject* pyFunction= PyDict_GetItemString(pDict, (char*)"energy");
     //PyObject* pyEnergy= PyObject_CallObject(pyFunction,NULL);
     PyObject* pyValue= PyObject_CallObject(pyFunction,NULL);
     double energy, density;

     PyArg_ParseTuple(pyValue, "dd", &energy, &density);
     //double energy= PyFloat_AsDouble(pyEnergy);

     Py_DECREF(pyModule);
     Py_DECREF(pyModuleString);

     return std::make_pair(density,energy);
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
     for(const auto& [b,s,include] : bs)
         bmax= max(bmax,b.cartesian().norm());

     int numberOfIgnoredPoints= 0;
     for (const auto &latticeVector: config) {
         VectorDimD x;
         OrderedTuplet<dim+1> temp;
         if (&(latticeVector.lattice) == &(this->gb.bc.A)) {
             double height= latticeVector.cartesian().dot(gb.nA.cartesian().normalized());
             if (height<FLT_EPSILON)
                 temp <<  gb.bc.getLatticeVectorInD(latticeVector),1;
             else
                 temp <<  gb.bc.getLatticeVectorInD(latticeVector),-1;
         }
         else if (&(latticeVector.lattice) == &(this->gb.bc.B)) {
             double height = latticeVector.cartesian().dot(gb.nB.cartesian().normalized());
             if (height<FLT_EPSILON)
                 temp <<  gb.bc.getLatticeVectorInD(latticeVector),2;
             else
                 temp <<  gb.bc.getLatticeVectorInD(latticeVector),-2;
         }

         //x = latticeVector.cartesian() + this->displacement(latticeVector.cartesian());

         if (&(latticeVector.lattice) == &(this->gb.bc.A)) {
             x = latticeVector.cartesian() + this->displacement(temp) + this->uAverage;
         }
         else if (&(latticeVector.lattice) == &(this->gb.bc.B) )
             x = latticeVector.cartesian() + this->displacement(temp) - this->uAverage;
         else
             x = latticeVector.cartesian() + this->displacement(temp);

         // ignore x if it occupies a deleted CSL position
         bool ignore= false;
         VectorDimD cslShift;
         cslShift << -0.5, -FLT_EPSILON, -FLT_EPSILON;
         VectorDimD xModulo= x;
         std::vector<LatticeVector<3>> localBoxVectors(boxVectors);
         localBoxVectors[0]=5*boxVectors[0];
         LatticeVector<dim>::modulo(xModulo, localBoxVectors, cslShift);
         for(const auto& [b,s, include] : bs) {
             if (include == 2 && (s - xModulo).norm() < 1e-6) {
                 ignore = true;
                 numberOfIgnoredPoints++;
                 break;
             }
         }
         if (ignore==true) {
             continue;
         }


         if (&(latticeVector.lattice) == &(this->gb.bc.A) && x.dot(this->gb.nA.cartesian().normalized()) <= 1e-6)
         //if (&(latticeVector.lattice) == &(this->gb.bc.A) && x.dot(this->gb.nA.cartesian().normalized()) <= FLT_EPSILON)
         //if (&(latticeVector.lattice) == &(this->gb.bc.A))
         {
             referenceConfigA.push_back(latticeVector.cartesian());
             deformedConfigA.push_back(x);
         }
         else if (&(latticeVector.lattice) == &(this->gb.bc.B) && x.dot(this->gb.nB.cartesian().normalized()) <= 1e-6)
         //else if (&(latticeVector.lattice) == &(this->gb.bc.B) && x.dot(this->gb.nB.cartesian().normalized()) <= FLT_EPSILON)
         //else if (&(latticeVector.lattice) == &(this->gb.bc.B))
         {
             referenceConfigB.push_back(latticeVector.cartesian());
             deformedConfigB.push_back(x);
         }

     }


     int numberOfIgnoredCSLPoints= 0;
     for(const auto& [b,s, include] : bs)
     {
         if (include==2) numberOfIgnoredCSLPoints++;
     }
     if(numberOfIgnoredPoints != 2*numberOfIgnoredCSLPoints)
         throw(std::runtime_error("GB Mesostate construction failed: incorrect number of points"));


     int nAtoms= referenceConfigA.size()+referenceConfigB.size()+configDscl.size();

     std::string referenceFile= name + "_reference0.txt";
     std::string deformedFile= name + "_reference1.txt";

     std::ofstream reference, deformed;
     reference.open(referenceFile);
     deformed.open(deformedFile);
     if (!reference || !deformed) std::cerr << "Unable to open files";
     reference << nAtoms << std::endl; deformed << nAtoms << std::endl;
     reference << "Lattice=\" "; deformed << "Lattice=\" ";

     reference << std::setprecision(15) << (2*boxVectors[0].cartesian()).transpose() << " ";
     deformed << std::setprecision(15) << (2*boxVectors[0].cartesian()).transpose() << " ";
     reference << std::setprecision(15) << (boxVectors[1].cartesian()).transpose() << " ";
     deformed << std::setprecision(15) << (boxVectors[1].cartesian()).transpose() << " ";
     reference << std::setprecision(15) << (boxVectors[2].cartesian()).transpose();
     deformed << std::setprecision(15) << (boxVectors[2].cartesian()).transpose();
     reference << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\" "; deformed << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\" F T T\" origin=\" ";
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
