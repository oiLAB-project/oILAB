//
// Created by Nikhil Chandra Admal on 2/4/24.
//
#include <MesoState.h>
#include <randomInteger.h>

namespace gbLAB
{
    template<int dim>
    MesoState<dim>::MesoState(const Gb<dim>& gb,
                              const ReciprocalLatticeVector<dim>& axis,
                              const int& periodScaling,
                              const double& a2,
                              const int& nImages):
            ReferenceState<dim>(gb,axis,periodScaling),
            Dislocations(a2,periodScaling*gb.getPeriodVector(axis).cartesian().norm(),nImages),
            currentState(this->refState)
    {
        /*
        const auto temp= getOrthogonalPlaneIndices(1);
        numberOfInteractingPlanes= floor((double)numberOfPlanesOrthogonalToGB()/
                          *std::min_element(temp.begin(), temp.end()));
        assert(numberOfInteractingPlanes>0);
        localState.resize(numberOfInteractingPlanes);
        */

    }
    template<int dim>
    MesoState<dim>::MesoState(const ReferenceState<dim>& rS, const double& a2, const int& nImages):
            ReferenceState<dim>(rS),
            currentState(this->refState),
            Dislocations(a2,this->periodScaling*this->gb.getPeriodVector(this->axis).cartesian().norm(),nImages)
    {}

    template<int dim>
    void MesoState<dim>::removeDislocation(const Triplet& t)
    {
        for(auto& index : currentState)
        {
            if (index[0]>=t(0) && index[0]<=t(1))
                index[2] = index[2] + t(2) / 2;
        }

        // find index of Triplet in defectsIndies
        int index= std::find(defectsIndices.begin(), defectsIndices.end(),t)-defectsIndices.begin();

        // update defectsIndices and dislocation dipoles
        this->removeDislocationDipole(index);
        defectsIndices.erase(defectsIndices.begin() + index);
    }

    template<int dim>
    Triplet MesoState<dim>::removeRandomDislocation()
    {
        int index= random<int>(0,defectsIndices.size()-1);
        Triplet t= defectsIndices[index];
        removeDislocation(t);
        return t;
    }

    template<int dim>
    void MesoState<dim>::insertDislocation(const Triplet& dipole)
    {
        // update the current state
        for(auto& index : currentState)
        {
            if (index[0]>=dipole(0) && index[0]<=dipole(1))
                index[2] = index[2] - dipole(2) / 2;
        }

        // update defectsIndices
        defectsIndices.push_back(dipole);

        Matrix2d ends;
        ends.col(0)= Vector2d(this->shiftLatticeVector().cartesian().norm()*((double) dipole(0)-0.5)/this->numberOfPlanesOrthogonalToGB(),0);
        ends.col(1)= Vector2d(this->shiftLatticeVector().cartesian().norm()*((double) dipole(1)+0.5)/this->numberOfPlanesOrthogonalToGB(),0);

        Vector2d b;
        //Vector2d shiftVector2d;
        b << 0,dipole(2)*this->gb.nA.planeSpacing();
        //shiftVector2d << this->shiftLatticeVector().cartesian().norm(),0;
        insertDislocationDipole(ends,b);
    }

    template<int dim>
    std::set<int> MesoState<dim>::getOrthogonalPlaneIndices(const int& parallelPlaneIndex) const
    {
        std::set<int> orthogonalPlaneIndices;
        for(const auto& index: currentState) {
            if (index(2) == parallelPlaneIndex)
                orthogonalPlaneIndices.insert(index(0));
        }
        return orthogonalPlaneIndices;

    }

    template<int dim>
    Triplet MesoState<dim>::insertRandomDislocation(const int& dipoleSign)
    {
        std::set<int> cslPointsIndices(getOrthogonalPlaneIndices(0));
        // insert the CSL point corresponding to the end of the period vector
        cslPointsIndices.insert(this->numberOfPlanesOrthogonalToGB());
        std::set<int> firstPlanePointsIndices(getOrthogonalPlaneIndices(dipoleSign));
        if (firstPlanePointsIndices.empty())
            throw (std::runtime_error("Unable to insert dislocation dipole of sign " + std::to_string(dipoleSign)));

        // Look for CSL points neighboring firstPlanePointsIndices[0]
        int leftCslPointIndex, rightCslPointIndex;
        for(auto iter= cslPointsIndices.begin(); iter!= cslPointsIndices.end(); ++iter)
        {
            leftCslPointIndex= *iter;
            rightCslPointIndex= *next(iter,1);
            if (rightCslPointIndex > *firstPlanePointsIndices.begin())
                break;
        }

        // dipoleBegin - random integer in [0, firstPlanePointsIndices[0]]
        const int& dipoleBegin= random<int>(leftCslPointIndex+1,*firstPlanePointsIndices.begin());

        assert(rightCslPointIndex > *firstPlanePointsIndices.begin());
        const int& dipoleEnd= random<int>(*firstPlanePointsIndices.begin(),rightCslPointIndex-1);


        Triplet output(dipoleBegin,dipoleEnd,dipoleSign*2);
        insertDislocation(output);
        return output;
    }

    template<int dim>
    Triplet MesoState<dim>::insertRandomDislocation()
    {
        std::set<int> cslPointsIndices(getOrthogonalPlaneIndices(0));
        // insert the CSL point corresponding to the end of the period vector
        cslPointsIndices.insert(this->numberOfPlanesOrthogonalToGB());

        // Look for CSL points neighboring firstPlanePointsIndices[0]
        assert(*cslPointsIndices.begin() +1 < *prev(cslPointsIndices.end()));

        int dipoleCenter= 0;
        while(cslPointsIndices.find(dipoleCenter) != cslPointsIndices.end())
            dipoleCenter= random<int>(*cslPointsIndices.begin()+1,
                                          *prev(cslPointsIndices.end())-1);

        int leftCslPointIndex, rightCslPointIndex;
        for(auto iter= cslPointsIndices.begin(); iter!= prev(cslPointsIndices.end()); ++iter)
        {
            leftCslPointIndex= *iter;
            rightCslPointIndex= *next(iter,1);
            if (rightCslPointIndex > dipoleCenter)
                break;
        }

        // dipoleBegin - random integer in [0,dipoleCenter]
        const int& dipoleBegin= random<int>(leftCslPointIndex+1,dipoleCenter);

        assert(rightCslPointIndex > dipoleCenter);
        const int& dipoleEnd= random<int>(dipoleCenter,rightCslPointIndex-1);

        int dipoleSign= random<int>(-100,100) <= 0 ? -1 : 1;
        Triplet output(dipoleBegin,dipoleEnd,dipoleSign*2);
        insertDislocation(output);
        return output;
    }

    template<int dim>
    //template<int dm=dim>
    typename std::enable_if<dim==3,void>::type
    MesoState<dim>::box(const int& heightFactor,
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
        boxVectors.push_back(this->shiftLatticeVector());
        boxVectors.push_back(vectorAlongAxisC);
        assert(abs(boxVectors[1].cartesian().dot(boxVectors[2].cartesian())) < FLT_EPSILON);

        auto config= this->gb.bc.box(boxVectors,1,1);
        std::vector<VectorDimD> referenceConfigA, deformedConfigA;
        std::vector<VectorDimD> referenceConfigB, deformedConfigB;
        std::vector<VectorDimD> configDscl;

        for (const auto &latticeVector: config) {
            if (&(latticeVector.lattice) == &(this->gb.bc.A) || &(latticeVector.lattice) == &(this->gb.bc.B))
            {
                VectorDimD coordinates;
                coordinates(0) = latticeVector.cartesian().dot(this->shiftLatticeVector().cartesian().normalized());
                coordinates(1) = latticeVector.cartesian().dot(this->gb.nA.cartesian().normalized());
                coordinates(2) = latticeVector.cartesian().dot(this->axis.cartesian().normalized());
                auto X2d(Eigen::Vector2d(coordinates(0), coordinates(1)));
                Vector2d x2d;

                try {
                    if ((&(latticeVector.lattice) == &(this->gb.bc.A) && latticeVector.dot(this->gb.nA) < 0) ||
                        (&(latticeVector.lattice) == &(this->gb.bc.B) && latticeVector.dot(this->gb.nB) < 0))
                        x2d = deformationMap(X2d, 1);
                    else
                        x2d = deformationMap(X2d, -1);
                }
                catch (std::runtime_error& e){
                    std::cout << e.what() << std::endl;
                    throw (std::runtime_error("Unable to calculate the deformation map"));
                }


                VectorDimD x3d;
                x3d = x2d(0) * this->shiftLatticeVector().cartesian().normalized() +
                      x2d(1) * this->gb.nA.cartesian().normalized() +
                      coordinates(2) * boxVectors[2].cartesian().normalized();

                if (&(latticeVector.lattice) == &(this->gb.bc.A) && x3d.dot(this->gb.nA.cartesian().normalized()) <= 1e-5)
                    //if (&(latticeVector.lattice) == &(gb.bc.A))
                {
                    referenceConfigA.push_back(latticeVector.cartesian());
                    deformedConfigA.push_back(x3d);
                }
                else if (&(latticeVector.lattice) == &(this->gb.bc.B) && x3d.dot(this->gb.nB.cartesian().normalized()) <= 1e-5)
                    //else if (&(latticeVector.lattice) == &(gb.bc.B))
                {
                    referenceConfigB.push_back(latticeVector.cartesian());
                    deformedConfigB.push_back(x3d);
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



    template<int dim>
    bool MesoState<dim>::operator<(const MesoState& rhs) const
    {
        assert(&(this->gb) == &(rhs.gb) &&
               &(this->axis)==&(rhs.axis) &&
               this->currentState.size() == rhs.currentState.size()) ;

        //if (numberOfCslPointsInGb < rhs.numberOfCslPointsInGb)
        if (getOrthogonalPlaneIndices(0).size()< rhs.getOrthogonalPlaneIndices(0).size())
            return true;
        //else if (rhs.numberOfCslPointsInGb < numberOfCslPointsInGb)
        else if (rhs.getOrthogonalPlaneIndices(0).size() < getOrthogonalPlaneIndices(0).size())
            return false;
        else // both have equal number of dislocations
        {
            for (int index = 0; index < rhs.currentState.size(); ++index) {
                if (!(this->currentState[index] < rhs.currentState[index]) &&
                    !(rhs.currentState[index] < this->currentState[index]))
                    continue;
                else if (this->currentState[index] < rhs.currentState[index])
                    return true;
                else
                    return false;
            }
        }
        return false;

    }


    template<int dim>
    double MesoState<dim>::energy() const
    {
        return elasticEnergy() + this->planeEnergies.dot(getLocalStateCount(10).template cast<double>());
    }


    template<int dim>
    template <typename T>
    int MesoState<dim>::sgn(T val) const
    {
        return (T(0) < val) - (val < T(0));
    }

    template class MesoState<3>;

}
