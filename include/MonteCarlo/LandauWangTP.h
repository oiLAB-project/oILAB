//
// Created by Nikhil Chandra Admal on 8/14/24.
//

#ifndef OILAB_LANDAUWANGTP_H
#define OILAB_LANDAUWANGTP_H
#include<EvolutionAlgorithm.h>
#include<vector>
#include<Eigen/Eigen>
#include <fstream>

namespace gbLAB {

// LandauWangTP is an EvolutionAlgorithm with an evolving transition probability
// include a switch which indicates whether the transition probability is allowed to evolve
    template<typename StateType, typename SystemType>
    class LandauWangTP : public EvolutionAlgorithm<StateType, SystemType, LandauWangTP<StateType,SystemType>> {
    private:
        bool exponentialRegime;
        double f;
        int countLW;
        double currentEnergy, currentDensity;
        const std::tuple<double,double,int> energyLimits, densityLimits;
        const int numberOfEnergyStates, numberOfDensityStates;
        Eigen::MatrixXi histogram;
        std::map<StateType, std::pair<double,double>> stateDensityEnergyMap;
        std::ofstream spectrumFile;


        bool histogramIsFlat(const double& c) const;

        static std::tuple<int,int,bool> spectrumIndex(const double& energy,
                                                      const double& density,
                                                      const std::tuple<double,double,int>& energyLimits,
                                                      const std::tuple<double,double,int>& densityLimits);
        static Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> getMask(const int& numberOfEnergyStates,
                                                                         const int& numberOfDensityStates);
        static std::map<StateType,std::pair<double,double>> getStateDensityEnergyMap();
        static Eigen::MatrixXd getTheta(const Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>& mask, double& f);

    public:
        Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> mask;
        Eigen::MatrixXd theta;


        explicit LandauWangTP(const std::tuple<double,double,int>& energyLimits);

        LandauWangTP(const std::tuple<double,double,int>& energyLimits,
                     const std::tuple<double,double,int>& densityLimits);

        double probability(const std::pair<StateType,SystemType>& proposedState,
                           const std::pair<StateType,SystemType>& currentState);

        void writeTheta(const std::string& filename) const;


    };
}
#endif //OILAB_LANDAUWANGTP_H
