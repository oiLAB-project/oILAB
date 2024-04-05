//
// Created by Nikhil Chandra Admal on 3/25/24.
//
#include <MesoStateEnsemble.h>
#include <randomInteger.h>
namespace gbLAB {
    template<int dim>
    MesoStateEnsemble<dim>::MesoStateEnsemble(const ReferenceState<dim>& rS, const std::string& filename):
    rS(rS)
    {
        std::ifstream mesostatesInFile;
        int numberOfDipoles= 0;
        std::string line;
        mesostatesInFile.open(filename);
        if ( mesostatesInFile.is_open() ) {
            while (mesostatesInFile.peek() != EOF ) {
                std::getline(mesostatesInFile, line);
                if (line.substr(0, line.find(" ")) == "Mesostate:")
                    continue;
                if (line.substr(0, line.find(" ")) == "Number") {
                    size_t lastIndex = line.find_last_not_of("0123456789");
                    numberOfDipoles= std::stoi(line.substr(lastIndex + 1));
                }
                MesoState<3> ms(rS,1,100);
                for(int i=0; i<numberOfDipoles; ++i)
                {
                    Triplet t;
                    for(int j=0;j<3;++j) {
                        if (!(mesostatesInFile >> t(j))) std::cout << "error reading" << std::endl;
                    }
                    ms.insertDislocation(t);
                }
                //this->insert(ms);
                this->push_back(ms);
                std::getline(mesostatesInFile, line);
            }
        }
        mesostatesInFile.close();
    }

    template<int dim>
    MesoStateEnsemble<dim>::MesoStateEnsemble(const ReferenceState<dim>& rS,
                                              const int& maxNumberOfMesoStates): rS(rS)
    {

        const int maxIterations= 10*maxNumberOfMesoStates;
        const int maxNumberOfDipoles= 10;
        int numberOfMesoStates=0;
        for(int i=0; i< maxIterations; ++i) {
            MesoState<3> ms(rS,1,100);
            for(int j=0; j< maxNumberOfDipoles && numberOfMesoStates < maxNumberOfMesoStates; ++j) {
                int dipoleSign= random<int>(-100,100) < 0 ? -1 : 1;
                try {
                    ms.insertRandomDislocation(dipoleSign);
                    numberOfMesoStates++;
                    //this->insert(ms);
                    this->push_back(ms);
                }
                catch(std::runtime_error& e) {
                    dipoleSign=-1*dipoleSign;
                    try {
                        ms.insertRandomDislocation(dipoleSign);
                        numberOfMesoStates++;
                        //this->insert(ms);
                        this->push_back(ms);
                    }
                    catch(std::runtime_error& e) {
                        break;
                    }
                }
            }
        }

    }

    template<int dim>
    MesoStateEnsemble<dim>::MesoStateEnsemble(const ReferenceState<dim>& rS): rS(rS)
    {
        MesoState<3> ms(rS,1,100);
        this->push_back(ms);
    }

    template<int dim>
    MesoStateEnsemble<dim> MesoStateEnsemble<dim>::build() const
    {
        MesoStateEnsemble<dim> builtMesostates(rS);
        int index= 0;
        int parentIndex= 0;
        for(const auto& ms : *this)
        {
            try {
                ms.box(2,1, "ms" + std::to_string(index));
                //builtMesostates.insert(ms);
                builtMesostates.push_back(ms);
                index++;
            }
            catch(std::runtime_error& e)
            {
                std::cout << e.what() << std::endl;
                std::cout << "Skipping mesostate " << parentIndex << std::endl;
                std::cout <<  std::endl;
            }
            parentIndex++;
        }
        return builtMesostates;
    }

    template<int dim>
    void MesoStateEnsemble<dim>::write(const std::string& filename) const
    {
        std::ofstream mesostatesFile;
        mesostatesFile.open(filename);
        int mesoStateIndex= 0;
        for(const auto& ms : *this) {
            mesostatesFile << "Mesostate: " << mesoStateIndex << std::endl;
            mesostatesFile << "Number of dipoles =  " << ms.defectsIndices.size() << std::endl;
            for (const auto& dipoleIndices: ms.defectsIndices)
                mesostatesFile << dipoleIndices.transpose() << std::endl;
            mesoStateIndex++;
        }
        mesostatesFile.close();
    }

    template class MesoStateEnsemble<3>;

}
