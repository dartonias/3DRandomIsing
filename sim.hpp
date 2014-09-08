#ifndef SIMHPP
#define SIMHPP

#ifndef DEBUG
#define DEBUG 0
#endif //DEBUG
// Stephen Inglis, 2014.06.23
// sim.hpp, header file for the main simulation part of the Random 3D Ising model
// Will eventually have two replicas so we can calculate the entanglement entropy

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "MersenneTwister.h"
#include "Eigen/Core"
#include "stat.hpp"
#include <assert.h>

#define PI 3.1415926535897

// Some helper functions, for flow control

bool fexists(const char* filename){
    ifstream ifile(filename);
    return ifile;
}

void getTemps(vector<double>& temps){
    std::string filename = "temps.dat";
    std::fstream inFile(filename.c_str());
    temps.resize(0);
    double cur_temp;
    while (inFile >> cur_temp){
        temps.push_back(cur_temp);
    }
}

void flip(int& z){
    z = z*-1;
}

class Sim{
    private:
        static const int binSize = 1000; // Averages per line (larger --> smaller files)
        static const int bufferSize = 10; // Lines held per write (larger --> less often writes to file)
        Eigen::Matrix<int, Eigen::Dynamic, 2> cluster; // For the cluster update
        Eigen::Matrix<int, Eigen::Dynamic, 1> Jmat; // All elements +- 1
        Eigen::Matrix<int, Eigen::Dynamic, 2> spins; // All spins +- 1
        void loadParams(); // Loads parameters from param.dat
        MTRand* rand; // MersenneTwister pseudorandom number generator
        int seed; // Seed for the random number generator
        int L; // Size of the lattice
        int Nspins; // Derived from L above
        int Nbonds; // Derived from L above
        double P; // Probability of disorder bond, 0 = pure ferromagnet
        double beta; // Inverse temperature
        int regionA; // Number of half sheets in regionA
        double Padd; // Probability of adding to a cluster
        double Energy; // Total energy
        int reverse; // Variable for calculating the regions in reverse size, since we lack symmetry
        int getAdj(int spin, int num); // Gets one of the 6 neighbours of a site
    public:
        Sim(double _beta=0.1); //Default constructor
        int sweeps; // Number of MC sweeps per bin
        int bins; // Number of total bins
        void singleUpdate(); // Single spin update attempt
        void saveJ();
        void loadJ();
        void addNeighbours(int z, int zz, int& s);
        void clusterUpdate(); // Swendsen-Wang update
};

int Sim::inA(int s){
    // Given a site, returns whether this site is in regionA, using the internally stored value of regionA
    // For each integer value of regionA, half a sheet of sites is added. So for an LxLxL system, regionA
    // varies from zero to 2L, where the number of layers in regionA is N/2
    if (reverse==1) s = L*L*L-1-s;
    if (s>((regionA/2+1)*L*L)){
        return 0;
    }
    else if(s<((regionA/2)*L*L)){
        return 1;
    }
    else if((regionA%2)==0){
        return 0;
    }
    // Interesting case, in the layer being added, and half the layer has been added so far
    // since regionA%2==1
    else if(((s/L)%2)==0){
        return 1;
    }
    else return 0;
}

Sim::Sim(double _beta){
    loadParams();
    beta = _beta;
    Padd = 1. - exp(-2*beta); // Probability of adding a spin when they are satisfied for the cluster move
    rand = new MTRand(seed);
    Nspins = L*L*L;
    Nbonds = 3*Nspins;
    spins.resize(Nspins,2);
    cluster.resize(Nspins,2)
    for(int i=0;i<Nspins;i++){
        spins(i,0) = rand->randInt(1)*2 - 1; // +- 1 random initial state
        if (inA(i)){
            spins(i,1) = spins(i,0); // spins must match in region A
        }
        else{
            spins(i,1) = rand->randInt(1)*2 - 1; // +- 1 random initial state
        }
    }
    Jmat.resize(Nbonds);
    for(int i=0;i<Nbonds;i++){
        if(rand->randExc() < P){
            Jmat(i) = 1; // If P is small, this usually won't happen
        }
        else{
            Jmat(i) = -1; // Usually we are ferromagnetic
        }
    }
    if(DEBUG){
        std::cout << Jmat.transpose() << std::endl;
    }
}

void Sim::loadParams(){
    std::string filename = "param.dat";
    std::string g; // Garbage string for going through param file
    std::fstream inFile(filename.c_str());
    inFile >> g >> L;
    inFile >> g >> P;
    inFile >> g >> seed;
    inFile >> g >> sweeps;
    inFile >> g >> bins;
    inFile >> g >> regionA;
    inFile >> g >> reverse;
    if(DEBUG){
        std::cout << "L = " << L << std::endl;
        std::cout << "P = " << P << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "seed = " << seed << std::endl;
        std::cout << "sweeps = " << sweeps << std::endl;
        std::cout << "bins = " << bins << std::endl;
        std::cout << "regionA = " << regionA << std::endl;
        std::cout << "reverse = " << reverse << std::endl;
    }
}

int Sim::getAdj(int s, int n){
    // Retrieves the n'th neighbour of spin s
    assert((n>=0)&&(n<6));
    assert((s>=0)&&(s<Nspins));
    if(n==0){
        // Spin right
        if((s%L)==(L-1)) return s+1-L;
        else return s+1;
    }
    else if(n==1){
        // Spin left
        if((s%L)==0) return s-1+L;
        else return s-1
    }
    else if(n==2){
        // Spin up
        if(((s/L)%L)==(L-1)) return s+L-L*L;
        else return s+L;
    }
    else if(n==3){
        // Spin down
        if(((s/L)%L)==0) return s-L+L*L;
        else return s-L;
    }
    else if(n==4){
        // Spin up one layer
        if(((s/L/L)%L)==(L-1)) return s+L*L-L*L*L;
        else return s+L*L;
    }
    else if(n==5){
        // Spin down one layer
        if(((s/L/L)%L)==0) return s-L*L+L*L*L;
        else return s-L*L;
    }
}

int Sim::getBond(int s, int n){
    // Retrieves the bond connecting the n'th neighbour of spin s
    assert((n>=0)&&(n<6));
    assert((s>=0)&&(s<Nspins));
    if(n==0){
        return 3*s+0;
    }
    else if(n==1){
        return 3*getAdj(s,1)+0;
    }
    else if(n==2){
        return 3*s+1;
    }
    else if(n==3){
        return 3*getAdj(s,3)+1;
    }
    else if(n==4){
        return 3*s+2;
    }
    else if(n==5){
        return 3*getAdj(s,5)+2;
    }
}

void Sim::singleUpdate(){
    // Update of a single spin
    int z = rand->randInt(Nspins-1); 
    double field = 0;
    if(inA(z)){
        for(int i=0;i<6;i++){
            field += spins(getAdj(z,i),0)*Jmat(getBond(z,i));
            field += spins(getAdj(z,i),1)*Jmat(getBond(z,i));
        }
        if(rand->randExc() < exp(-2*spins(z,0)*field)){
            flip(spins(z,0));
            flip(spins(z,1));
        }
    }
    else{
        int layer = rand->randInt(1); // Generates a 0 or 1
        for(int i=0;i<6;i++){
            field += spins(getAdj(z,i),layer)*Jmat(getBond(z,i));
        }
        if(rand->randExc() < exp(-2*spins(z,0)*field)){
            flip(spins(z,layer));
        }
    }
}

void Sim::addNeighbours(int z,int zz,int& s){ // Recursive part
    s = s+1;
    cluster(z,zz) = -1;
    if(inA(z)){
        if(cluster(z,(zz+1)%2)==1) addNeighbours(z,(zz+1)%2,s);
    }

    for(int i=0;i<6;i++){
        if(cluster(getAdj(z,i),zz)==1){ // Only try if it's not in the cluster
            if((spins(z,zz) * spins(getAdj(z,i),zz) * Jmat(getBond(z,i)))==-1){ // If the spins match ...
                if(rand->randExc() < Padd){ // Probability to add
                    addNeighbours(adjS(z,i),zz,s);
                }
            }
        }
    }
}

void Sim::clusterUpdate(){
    // Choose a random spin
    int z = rand->randInt(Nspins-1);
    int zz = rand->randInt(1);
    cluster.fill(1); // No spins are in the cluster
    if(DEBUG){
        int s = 0;
        addNeighbours(z,zz,s);
        std::cout << s << std::endl;
        for(int i=0;i<Nspins;i++){
            if((i%L)==0) std::cout << std::endl;
            std::cout << std::setw(2);
            std::cout << cluster(i) << " ";
        }
        std::cout << std::endl;
    }
    else{
        int s = 0;
        addNeighbours(z,zz,s);
    }
    // Now that it's finished, flip all spins in the cluster
    // All elements in "cluster" are -1, and will flip spins in "spins"
    if(DEBUG){
        updateE();
        double tE1 = Energy;
        spins = spins.array() * cluster.array();
        updateE();
        double tE2 = Energy;
        std::cout << "dE = " << (tE2 - tE1) << std::endl << std::endl;
    }
    else spins = spins.array() * cluster.array();
}

void Sim::saveJ(){
    ofstream jfile;
    jfile.open("Jmat.dat");
    for(int i=0;i<Nbonds;i++){
        jfile << Jmat(i) << " ";
    }
    jfile.close();
}

void Sim::loadJ(){
    ifstream jfile;
    jfile.open("Jmat.dat");
    for(int i=0;i<Nbonds;i++){
        jfile >> Jmat(i);
    }
    jfile.close();
}

void Sim::updateE(){
    Energy = 0;
    // Loop over spins
    for(int i=0;i<Nspins;i++){
        // Loop over unique bonds, 2 per spin
        for(int b=0;b<3;b++){
            Energy += spins(i,0) * spins(getAdj(i,b),0) * Jmat(getBond(i,b));
            Energy += spins(i,1) * spins(getAdj(i,b),1) * Jmat(getBond(i,b));
        }
    }
}

#endif //SIMHPP
