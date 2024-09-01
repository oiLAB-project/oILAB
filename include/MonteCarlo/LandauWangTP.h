//
// Created by Nikhil Chandra Admal on 8/14/24.
//

#ifndef OILAB_LANDAUWANGTP_H
#define OILAB_LANDAUWANGTP_H
#include<EvolutionAlgorithm.h>
#include<vector>
#include<Eigen/Eigen>

namespace gbLAB {

// LandauWangTP is an EvolutionAlgorithm with an evolving transition probability
// include a switch which indicates whether the transition probability is allowed to evolve
    template<typename StateType, typename SystemType>
    class LandauWangTP : public EvolutionAlgorithm<StateType, SystemType, LandauWangTP<StateType,SystemType>> {
    private:
        const int numberOfEnergyStates, numberOfDensityStates;
        const std::tuple<double,double,int> energyLimits, densityLimits;

        int countLW;
        bool exponentialRegime;
        double f;
        double currentEnergy, currentDensity;
        Eigen::MatrixXi histogram;

        bool histogramIsFlat(const double& c) const;
        static std::tuple<int,int,bool> spectrumIndex(const double& energy,
                                                      const double& density,
                                                      const std::tuple<double,double,int>& energyLimits,
                                                      const std::tuple<double,double,int>& densityLimits);
        static Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> getMask(const std::tuple<double,double,int>& energyLimits,
                                                                         const std::tuple<double,double,int>& densityLimits,
                                                                         const std::string& filename);

    public:
        Eigen::MatrixXd theta;
        std::map<StateType, std::pair<double,double>> stateDensityEnergyMap;
        Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> mask;


        LandauWangTP(const double& minEnergy,
                     const double& maxEnergy,
                     const int& numberOfEnergyStates,
                     const std::string& filename="");

        LandauWangTP(const double& minEnergy, const double& maxEnergy, const int& numberOfEnergyStates,
                     const double& minDensity, const double& maxDensity, const int& numberOfDensityStates,
                     const std::string& filename="");

        double probability(const std::pair<StateType,SystemType>& proposedState,
                           const std::pair<StateType,SystemType>& currentState);


    };
}
#endif //OILAB_LANDAUWANGTP_H
