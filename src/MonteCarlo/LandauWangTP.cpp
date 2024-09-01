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
/*
template<typename StateType>
double LandauWangTP<StateType>::f= exp(1.0);

template<typename StateType>
Eigen::VectorXd LandauWangTP<StateType>::theta;
template<typename StateType>
std::vector<int> LandauWangTP<StateType>::histogram;
template<typename StateType>
bool LandauWangTP<StateType>::exponentialRegime= true;
template<typename StateType>
int LandauWangTP<StateType>::countLW= 0;
 */
    template<typename StateType, typename SystemType>
    LandauWangTP<StateType,SystemType>::LandauWangTP(const double& minEnergy,
                                                     const double& maxEnergy,
                                                     const int& numberOfEnergyStates,
                                                     const std::string& filename):
            numberOfEnergyStates(numberOfEnergyStates),
            numberOfDensityStates(1),
            energyLimits(std::make_tuple(minEnergy,maxEnergy,numberOfEnergyStates)),
            densityLimits(std::make_tuple(0.0,1.0,1)),
            histogram(Eigen::MatrixXi::Zero(numberOfEnergyStates,numberOfDensityStates)),
            theta(Eigen::MatrixXd::Ones(numberOfEnergyStates,numberOfDensityStates)),
            mask(getMask(energyLimits,densityLimits,filename)),
            exponentialRegime(true),
            f(exp(1.0)),
            countLW(-1),
            currentEnergy(0.0),
            currentDensity(0.0)
    {
        double thetaNorm= 0.0;
        for(int i=0; i<histogram.rows(); ++i)
        {
            for(int j=0; j<histogram.cols(); ++j)
            {
                if(!mask(i,j))
                    thetaNorm+= abs(theta(i,j));
                else
                    theta(i,j)= 0.0;
            }

        }
        theta= theta/thetaNorm;
        std::cout << "Number of the mask-free histogram bins= " << histogram.size()-mask.template cast<int>().sum() << std::endl;
    }

    template<typename StateType, typename SystemType>
    LandauWangTP<StateType,SystemType>::LandauWangTP(const double& minEnergy, const double& maxEnergy, const int& numberOfEnergyStates,
                                                     const double& minDensity, const double& maxDensity, const int& numberOfDensityStates,
                                                     const std::string& filename):
            numberOfEnergyStates(numberOfEnergyStates),
            numberOfDensityStates(numberOfDensityStates),
            energyLimits(std::make_tuple(minEnergy,maxEnergy,numberOfEnergyStates)),
            densityLimits(std::make_tuple(minDensity,maxDensity,numberOfDensityStates)),
            histogram(Eigen::MatrixXi::Zero(numberOfEnergyStates,numberOfDensityStates)),
            theta(Eigen::MatrixXd::Ones(numberOfEnergyStates,numberOfDensityStates)),
            mask(getMask(energyLimits,densityLimits,filename)),
            exponentialRegime(true),
            f(exp(1.0)),
            countLW(-1),
            currentEnergy(0.0),
            currentDensity(0.0)
    {
        double thetaNorm= 0.0;
        for(int i=0; i<theta.rows(); ++i)
        {
            for(int j=0; j<theta.cols(); ++j)
            {
                if(!mask(i,j))
                    thetaNorm+= abs(theta(i,j));
                else
                    theta(i,j)= 0.0;
            }

        }
        theta= theta/thetaNorm;
        std::cout << "Number of the mask-free histogram bins= " << histogram.size()-mask.template cast<int>().sum() << std::endl;
    }
    template<typename StateType, typename SystemType>
    double LandauWangTP<StateType,SystemType>::probability(const std::pair<StateType,SystemType>& proposedStateSystem,
                                                           const std::pair<StateType,SystemType>& currentStateSystem)
    {
        countLW++;
        std::cout << countLW << ")" << std::endl;
        const auto& currentState= currentStateSystem.first;
        const auto& currentSystem= currentStateSystem.second;
        const auto& proposedState= proposedStateSystem.first;
        const auto& proposedSystem= proposedStateSystem.second;

        // Compute current state properties
        // if this is the first call, compute the current energy and density
        if(countLW==0) {
            const auto& temp= currentSystem.densityEnergy();
            currentDensity= temp.first;
            currentEnergy= temp.second;
            stateDensityEnergyMap[currentState]= temp;
        }
        else {
            assert(stateDensityEnergyMap.find(currentState) != stateDensityEnergyMap.end());
            currentDensity= stateDensityEnergyMap.at(currentState).first;
            currentEnergy= stateDensityEnergyMap.at(currentState).second;
        }

        // Compute proposed state properties
        double proposedEnergy, proposedDensity;
        if (stateDensityEnergyMap.find(proposedState) != stateDensityEnergyMap.end()) {
            proposedDensity= stateDensityEnergyMap.at(proposedState).first;
            proposedEnergy = stateDensityEnergyMap.at(proposedState).second;
        }
        else {
            const auto& temp= proposedSystem.densityEnergy();
            proposedDensity= temp.first;
            proposedEnergy= temp.second;
            stateDensityEnergyMap[proposedState]= temp;
            std::cout << "proposed density = " << proposedDensity << ", energy = " << proposedEnergy << std::endl;
        }

        auto [currentEnergyIndex,currentDensityIndex,currentInSpectrum]= spectrumIndex(currentEnergy,currentDensity,energyLimits,densityLimits);
        auto [proposedEnergyIndex,proposedDensityIndex,proposedInSpectrum]= spectrumIndex(proposedEnergy,proposedDensity,energyLimits,densityLimits);

        //if (currentEnergyIndex < 0) currentEnergyIndex = 0;

        // if current state is outside the energy-density spectrum (or masked) accept the move
        if (!currentInSpectrum || mask(currentEnergyIndex,currentDensityIndex))
            return 1.0;
        // now we are guaranteed that the current state is within the energy spectrum and is unmasked

        // if the proposed state is not in the density-energy spectrum or is masked reject it
        if (!proposedInSpectrum || mask(proposedEnergyIndex,proposedDensityIndex))
            return 0.0;

        // update histogram and theta only when the proposed state is unmasked and is in the spectrum
        histogram(currentEnergyIndex,currentDensityIndex) += 1;
        theta(currentEnergyIndex,currentDensityIndex) *= f;
        double thetaNorm = theta.lpNorm<1>();
        theta = theta / thetaNorm;

        // output histogram
        std::cout << "current density = " << currentDensity
                  << ", energy = " << currentEnergy << std::endl;
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
    Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> LandauWangTP<StateType,SystemType>::getMask(const std::tuple<double,double,int>& energyLimits,
                                                                                                  const std::tuple<double,double,int>& densityLimits,
                                                                                                  const std::string& filename)
    {
        const auto& [minEnergy,maxEnergy,numberOfEnergyStates]= energyLimits;
        const auto& [minDensity,maxDensity,numberOfDensityStates]= densityLimits;
        Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> output(numberOfEnergyStates,numberOfDensityStates);
        output.setZero();

        if (filename.empty())
            return output;
        else {
            output.setOnes();
            std::ifstream file;
            file.open(filename);
            if (!file) std::cerr << "Unable to open file";
            double energy, density;

            while (file >> density >> energy)
            {
                auto [energyIndex,densityIndex,inSpectrum]= spectrumIndex(energy,density,energyLimits,densityLimits);
                if(inSpectrum)
                    output(energyIndex,densityIndex)= false;
            }
        }

        std::cout << output << std::endl;
        return output;
    }

    template class LandauWangTP<XTuplet,GbMesoState<3>>;

}
