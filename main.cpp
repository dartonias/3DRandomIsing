#include "sim.hpp"
#include <vector>
#include <cmath>

int main(int argc, char** argv){
    Eigen::Matrix<int, Eigen::Dynamic, 2> temp_spins;
    double temp_E;
    vector<double> temps;
    getTemps(temps);
    int num_temps = temps.size();
    vector<Sim> simVec;
    simVec.resize(0);
    simVec.push_back(Sim(temps[0]));
    if(fexists("LOADJ")){
        simVec[0].loadJ();
    }
    if(fexists("SAVEJ")){
        simVec[0].saveJ();
        return 0;
    }
    for(int i=1;i<num_temps;i++){
        simVec.push_back(Sim(temps[i],simVec[0].getJ(),simVec[0].getRand()));
    }
    for(int i_bins=0;i_bins<simVec[0].bins;i_bins++){
        for(int i_sweeps=0;i_sweeps<simVec[0].sweeps;i_sweeps++){
            for(int z=0;z<num_temps;z++){
                for(int i_onesweep=0;i_onesweep<simVec[0].getNspins();i_onesweep++){
                    simVec[z].singleUpdate();
                }
                simVec[z].clusterUpdate();
                simVec[z].updateE(); // Must do this, needed for parallel tempering
                if (i_sweeps == 0){
                    simVec[z].resetFac();
                }
                if (i_sweeps == NUM_AVG){
                    simVec[z].printFac();
                }

                if (i_sweeps<NUM_AVG){
                    simVec[z].updateRatio(1);
                }
                else{
                    simVec[z].updateRatio();
                }
                simVec[z].updateBinder();
            }
            // Parallel tempering step
            //if(num_temps>1){
            if(false){
                for(int z=0;z<num_temps;z++){
                    int flipme = simVec[0].getRand()->randInt(num_temps-2);
                    double dE = simVec[flipme+1].getE() - simVec[flipme].getE();
                    double dB = simVec[flipme+1].getB() - simVec[flipme].getB();
                    if(simVec[0].getRand()->randExc() < exp(dE*dB)){
                        temp_spins = simVec[flipme].getSpins();
                        temp_E = simVec[flipme].getE();
                        simVec[flipme].setSpins(simVec[flipme+1].getSpins());
                        simVec[flipme].setE(simVec[flipme+1].getE());
                        simVec[flipme+1].setSpins(temp_spins);
                        simVec[flipme+1].setE(temp_E);
                    }
                }
            }
        }
    }
}
