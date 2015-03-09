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
    Sim tempSim;
    simVec.resize(0);
    simVec.push_back(Sim(temps[0]));
    #ifdef NEWT
    vector<double> simHist;
    simHist.assign(num_temps-1,0);
    vector<double> deltaB;
    deltaB.assign(num_temps-1,0);
    #endif
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
                #ifndef NEWT
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
                #endif
            }
            // Temperature modification step
            // Also record success rates to file, for debugging purposes
            #ifdef NEWT
            if(num_temps>1){
                for(int z=0;z<num_temps-1;z++){
                    double dE = simVec[z+1].getE() - simVec[z].getE();
                    double dB = simVec[z+1].getB() - simVec[z].getB();
                    double prob = exp(dE*dB);
                    if(prob > 1){
                        simHist[z] += 1;
                    }
                    else{
                        simHist[z] += prob;
                    }
                }
            }
            #endif
            // Parallel tempering step
            if(num_temps>1){
                for(int z=0;z<num_temps;z++){
                    int flipme = simVec[0].getRand()->randInt(num_temps-2);
                    double dE = simVec[flipme+1].getE() - simVec[flipme].getE();
                    double dB = simVec[flipme+1].getB() - simVec[flipme].getB();
                    if(simVec[0].getRand()->randExc() < exp(dE*dB)){
                        tempSim = simVec[flipme];
                        simVec[flipme] = simVec[flipme+1];
                        simVec[flipme+1] = tempSim;
                    }
                }
            }
        }
        // Modifying the temperatures according to the acceptance rates
        // Also, printing the average acceptance rate to file
        #ifdef NEWT
        ofstream newTemps;
        newTemps.open("temps2.dat");
        ofstream acceptance;
        acceptance.open("acceptance.dat",std::ofstream::app);
        double sum = 0;
        if(num_temps>1){
            for(int z=0;z<num_temps-1;z++){
                simHist[z] /= simVec[0].sweeps;
                // Prevent too much change per step
                simHist[z] = ((simHist[z] > 0.1) ? simHist[z] : 0.1);
                // These numbers are all smaller than 1 before taking the square root, so they all become closer to 1
                simHist[z] = pow(simHist[z],0.5);
                sum += simHist[z];
                deltaB[z] = simVec[z+1].getB() - simVec[z].getB();
            }
            sum /= (num_temps-1);
            // sum is now strictly less than 1, as before it was strictly less than num_temps-1
            double fulldB = simVec[num_temps-1].getB() - simVec[0].getB();
            double tdB = 0;
            for(int z=0;z<num_temps-1;z++){
                simHist[z] /= sum;
                tdB += deltaB[z] * simHist[z];
            }
            for(int z=0;z<num_temps-1;z++){
                simHist[z] /= (tdB / fulldB);
                simVec[z+1].setB(simVec[z].getB() + deltaB[z]*simHist[z]);
                newTemps << simVec[z].getB() << std::endl;
            }
            newTemps << simVec[num_temps-1].getB() << std::endl;
            std::cout << sum << std::endl;
            acceptance << sum << std::endl;
        }
        simHist.assign(num_temps-1,0);
        newTemps.close();
        acceptance.close();
        #endif
    }
}
