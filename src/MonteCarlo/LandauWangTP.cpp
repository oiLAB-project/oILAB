//
// Created by Nikhil Chandra Admal on 8/14/24.
//
#include <LandauWangTP.h>
#include <iostream>
#include <numeric>
#include <OrderedTuplet.h>
#include <GbMesoState.h>
#include <iomanip>

namespace gbLAB {
/* ---------------------------------------------------*/
    template<typename StateType, typename SystemType>
    LandauWangTP<StateType,SystemType>::LandauWangTP(const std::tuple<double,double,int>& energyLimits):
            LandauWangTP(energyLimits,{0.0,1.0,1})
    {}

    template<typename StateType, typename SystemType>
    LandauWangTP<StateType,SystemType>::LandauWangTP(const std::tuple<double,double,int>& energyLimits,
                                                     const std::tuple<double,double,int>& densityLimits) try:
            exponentialRegime(true),
            f(exp(1.0)),
            countLW(-1),
            currentEnergy(0.0),
            currentDensity(0.0),
            energyLimits(energyLimits),
            densityLimits(densityLimits),
            numberOfEnergyStates(std::get<2>(energyLimits)),
            numberOfDensityStates(std::get<2>(densityLimits)),
            histogram(Eigen::MatrixXi::Zero(numberOfEnergyStates,numberOfDensityStates)),
            stateDensityEnergyMap(getStateDensityEnergyMap()),
            mask(getMask(numberOfEnergyStates,numberOfDensityStates)),
            theta(getTheta(mask,f))
    {
        std::cout << "Number of the mask-free histogram bins= " << histogram.size()-mask.template cast<int>().sum() << std::endl;
        spectrumFile.open("energyDensityLW.txt",std::ios_base::app);
    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        exit(0);
    }
    template<typename StateType, typename SystemType>
    double LandauWangTP<StateType,SystemType>::probability(const std::pair<StateType,SystemType>& proposedStateSystem,
                                                           const std::pair<StateType,SystemType>& currentStateSystem)
    {
        countLW++;
        std::cout << countLW << ") ";
        const auto& currentState= currentStateSystem.first;
        const auto& currentSystem= currentStateSystem.second;
        const auto& proposedState= proposedStateSystem.first;
        const auto& proposedSystem= proposedStateSystem.second;

        // Compute current state properties
        // if this is the first call, compute the current energy and density
        if(countLW==0) {
            const auto& temp= currentSystem.densityEnergy();
            //currentDensity= temp.first;
            currentDensity= currentState.density();
            currentEnergy= temp.second;
            std::cout << currentState;
            stateDensityEnergyMap[currentState]= temp;
            spectrumFile << currentDensity << " " << currentEnergy << " " << currentState << std::endl;
        }
        else {
            assert(stateDensityEnergyMap.find(currentState) != stateDensityEnergyMap.end());
            currentDensity= stateDensityEnergyMap.at(currentState).first;
            currentEnergy= stateDensityEnergyMap.at(currentState).second;
        }

        // Compute proposed state properties
        double proposedEnergy, proposedDensity;
        if (stateDensityEnergyMap.find(proposedState) != stateDensityEnergyMap.end()) {
            std::cout << "found. Proposed state = " << proposedState << std::endl;
            proposedDensity= stateDensityEnergyMap.at(proposedState).first;
            proposedEnergy = stateDensityEnergyMap.at(proposedState).second;
            std::cout << "proposedDensity =  " << proposedDensity << "; proposedEnergy = " << proposedEnergy << std::endl;
        }
        else {
            std::cout << "new" << std::endl;
            const auto& temp= proposedSystem.densityEnergy();
            //proposedDensity= temp.first;
            proposedDensity= proposedState.density();
            proposedEnergy= temp.second;
            stateDensityEnergyMap[proposedState]= temp;
            spectrumFile << proposedDensity << " " << proposedEnergy << " " << proposedState << std::endl;
        }

        auto [currentEnergyIndex,currentDensityIndex,currentInSpectrum]= spectrumIndex(currentEnergy,currentDensity,energyLimits,densityLimits);
        auto [proposedEnergyIndex,proposedDensityIndex,proposedInSpectrum]= spectrumIndex(proposedEnergy,proposedDensity,energyLimits,densityLimits);

        //if (currentEnergyIndex < 0) currentEnergyIndex = 0;

        // if current state is outside the energy-density spectrum (or masked) accept the move
        if (!currentInSpectrum || mask(currentEnergyIndex,currentDensityIndex))
            return 1.0;
        // now we are guaranteed that the current state is within the energy spectrum and is unmasked

        // if the proposed state is not in the density-energy spectrum or is masked reject it
        if (!proposedInSpectrum || mask(proposedEnergyIndex, proposedDensityIndex))
            return 0.0;

        // update histogram and theta only when the proposed state is unmasked and is in the spectrum
        histogram(currentEnergyIndex,currentDensityIndex) += 1;
        theta(currentEnergyIndex,currentDensityIndex) *= f;
        double thetaNorm = theta.lpNorm<1>();
        theta = theta / thetaNorm;

        // output histogram
        //std::cout << "current density = " << currentDensity
        //          << ", energy = " << currentEnergy << std::endl;
        for(int i=0; i<histogram.rows(); ++i)
        {
            for(int j=0; j<histogram.cols(); ++j)
            {
                if(!mask(i,j))
                    std::cout << histogram(i,j) << " ";
            }
        }
        std::cout << std::endl;

        // update f
        if (exponentialRegime) {
            if (histogramIsFlat(0.8)) {
                std::cout << "Histogram is flat: f = " << f << std::endl;
                // reset the histogram
                histogram.setZero();
                // exponential regime
                f = sqrt(f);
            }
            if (f < exp(1.0 / (countLW + 1))) {
                std::cout << "Beginning non-exponential regime." << std::endl;
                exponentialRegime = false;
                // non-exponential regime
                f = exp(1.0 / (countLW + 1));
            }
        } else
            f = exp(1.0 / (countLW + 1));

        return theta(currentEnergyIndex,currentDensityIndex) / theta(proposedEnergyIndex,proposedDensityIndex);
    }


    template<typename StateType,typename SystemType>
    bool LandauWangTP<StateType,SystemType>::histogramIsFlat(const double &c) const
    {
        //int sum = std::accumulate(histogram.begin(), histogram.end(), 0);
        int sum = histogram.sum();
        int size= histogram.size()-mask.template cast<int>().sum();

        std::vector<double> flatnessMeasure;

        /*
        for (const auto& element : histogram.reshaped())
            flatnessMeasure.push_back(abs((double) element / sum - 1. / size));
        // max(|e-mean|)/sum * size/0.8 < 1 => max(|e-mean|) < 0.8 * mean => max(mean-e) < 0.8 *mean => mean-min(e) < 0.8*mean => min(e) > 0.2*mean
         */

        for(int i=0; i<histogram.rows(); ++i)
        {
            for(int j=0; j<histogram.cols(); ++j)
            {
                if(!mask(i,j))
                    flatnessMeasure.push_back(abs((double) histogram(i,j)/ sum - 1. / size));
            }

        }

        double nonFlatness= *std::max_element(flatnessMeasure.begin(), flatnessMeasure.end())*size/c;
        //if (*std::max_element(flatnessMeasure.begin(), flatnessMeasure.end()) < c / histogram.size())
        std::cout << "Flatness measure = " << nonFlatness << std::endl;

        if (nonFlatness<1)
            return true;
        else
            return false;
    }

    template<typename StateType,typename SystemType>
    std::tuple<int,int,bool> LandauWangTP<StateType,SystemType>::spectrumIndex(const double& energy,
                                                                               const double& density,
                                                                               const std::tuple<double,double,int>& energyLimits,
                                                                               const std::tuple<double,double,int>& densityLimits)
    {
        const auto& [minEnergy,maxEnergy,numberOfEnergyStates]= energyLimits;
        const auto& [minDensity,maxDensity,numberOfDensityStates]= densityLimits;
        double deltaE= (maxEnergy - minEnergy) / numberOfEnergyStates;
        double deltarho= (maxDensity - minDensity) / numberOfDensityStates;

        int energyIndex = floor((energy - minEnergy) / deltaE);
        int densityIndex = floor((density- minDensity) / deltarho);
        if (energyIndex >= 0 && energyIndex < numberOfEnergyStates &&
            densityIndex >= 0 && densityIndex < numberOfDensityStates)
            return std::make_tuple(energyIndex,densityIndex,true);
        else
            return std::make_tuple(energyIndex,densityIndex,false);
    }

    template<typename StateType,typename SystemType>
    Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> LandauWangTP<StateType,SystemType>::getMask(const int& numberOfEnergyStates,
                                                                                                  const int& numberOfDensityStates)
    {
        Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> output(numberOfEnergyStates,numberOfDensityStates);
        output.setZero();
        std::ifstream file;
        file.open("mask.txt");
        std::string line;
        if (!file)
            std::cout << "File mask.txt not found.";
        else{
            int i= 0;
            while(std::getline(file,line)) {
                std::stringstream s(line);
                for (int j = 0; j < numberOfDensityStates; ++j)
                    s >> output(i, j);
                ++i;
            }
        }
        std::cout << "mask = " << std::endl;
        std::cout << output << std::endl;
        return output;
    }

    template<typename StateType,typename SystemType>
    std::map<StateType, std::pair<double,double>> LandauWangTP<StateType,SystemType>::getStateDensityEnergyMap()
    {
        std::map<StateType, std::pair<double,double>> output;
        std::ifstream file;
        file.open("energyDensityLW.txt");
        if (!file)
            return output;

        double energy, density;

        std::string line;
        while(std::getline(file,line)) {
            std::stringstream s(line);
            s >> density; s >> energy;
            int temp;
            std::vector<int> tempVector;
            while (s >> temp)
                tempVector.push_back(temp);

            StateType state(tempVector.size());
            for (int i=0; i<tempVector.size(); ++i)
                state(i)= tempVector[i];

            output[state]= std::make_pair(density,energy);
        }

        /*
        for (const auto& [key,value]: output)
            std::cout <<  value.first << " " << value.second << " " << key << std::endl;
         */
        std::cout << "Number of pre-computed states = " << output.size() << std::endl;
        file.close();
        return output;
    }

    template<typename StateType,typename SystemType>
    Eigen::MatrixXd LandauWangTP<StateType,SystemType>::getTheta(const Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>& mask,
                                                                 double& f)
    {
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> outputTheta(mask.rows(),mask.cols());
        std::ifstream file;
        file.open("theta.txt");
        std::string line;
        if (!file) {
            std::cout << "File theta.txt not found. Setting a uniform distribution for theta";
            outputTheta.setOnes();
        }
        else {
            outputTheta.setZero();
            int lno= -1;
            while(std::getline(file,line)) {
                std::stringstream s(line);
                int j= 0;
                double temp;
                while(s >> temp)
                {
                    if (lno==-1) {
                        f = temp;
                        break;
                    }
                    outputTheta(lno,j)= temp;
                    j++;
                }
                if (lno != -1 && j != mask.cols())
                    throw(std::runtime_error("Error reading row in theta.txt") );
                lno++;
            }
            if (lno != mask.rows())
                throw(std::runtime_error("Error reading column in theta.txt") );

        }

        // normalize theta
        double outputThetaNorm= 0.0;
        for(int i=0; i<outputTheta.rows(); ++i)
        {
            for(int j=0; j<outputTheta.cols(); ++j)
            {
                if(!mask(i,j))
                    outputThetaNorm+= abs(outputTheta(i,j));
                else
                    outputTheta(i,j)= 0.0;
            }
        }
        outputTheta= outputTheta/outputThetaNorm;
        std::cout << "f = " << f << std::endl;
        std::cout << "theta = " << std::endl;
        std::cout << outputTheta << std::endl;
        return outputTheta;
    }


    template<typename StateType,typename SystemType>
    void LandauWangTP<StateType,SystemType>::writeTheta(const std::string& filename) const
    {
        std::ofstream outputFileHandle;
        outputFileHandle.open(filename);
        outputFileHandle << f << std::endl;
        outputFileHandle << theta;
        outputFileHandle.close();
    }

    template class LandauWangTP<XTuplet,GbMesoState<3>>;

}
