#ifndef DMRG_H
#define DMRG_H

#include "itensor/all.h"
#include <cstdlib>
#include <random>
#include <fstream>

using namespace std;
using namespace itensor;

class DMRG {
   
    // TODO: need an MPS collapse function
    int N_;
    MPO Hamiltonian_;
    SiteSet sites_;
    MPS psi_;
    MPS psi0_;

    public:

    DMRG(int N) : N_(N) {}

    inline MPS GetPsi() {
        return psi_;
    }

    inline SiteSet GetSiteSet() {
        return sites_;
    }

    void LTFIM1D(double J, double h, double Omega, string BC) {
        sites_ = SpinHalf(N_, {"ConserveQNs=",false});
        auto ampo = AutoMPO(sites_);
       
        for(int j = 1; j < N_; ++j) {
            ampo += -J*4.0,"Sz",j,"Sz",j+1;
            ampo += -h*2.0,"Sx",j;
            ampo += -Omega*2.0,"Sz",j;
        }
        ampo += -h*2.0,"Sx",N_;
        ampo += -Omega*2.0,"Sz",N_;
    
        if (BC == "PBC") {
            ampo += -J*4.0,"Sz",1,"Sz",N_;
        }
    
        Hamiltonian_ = toMPO(ampo);
    }
    
    void RotTFIM1D(double J, double h, string BC) {

        sites_ = SpinHalf(N_);
        auto ampo = AutoMPO(sites_);
       
        for(int j = 1; j < N_; ++j) {
            ampo += -J*4.0,"Sx",j,"Sx",j+1;
            ampo += -h*2.0,"Sz",j;
        }
        ampo += -h*2.0,"Sz",N_;
    
        if (BC == "PBC") {
            ampo += -J*4.0,"Sx",1,"Sx",N_;
        }
    
        Hamiltonian_ = toMPO(ampo);
    }
    
    void TFIM1D(double J, double h, string BC) {

        sites_ = SpinHalf(N_, {"ConserveQNs=",false});
        auto ampo = AutoMPO(sites_);
       
        for(int j = 1; j < N_; ++j) {
            ampo += -J*4.0,"Sz",j,"Sz",j+1;
            ampo += -h*2.0,"Sx",j;
        }
        ampo += -h*2.0,"Sx",N_;
    
        if (BC == "PBC") {
            ampo += -J*4.0,"Sz",1,"Sz",N_;
        }
    
        Hamiltonian_ = toMPO(ampo);
    }
    
    void Heisenberg1D(double J, string BC) {

        sites_ = SpinHalf(N_);
        auto ampo = AutoMPO(sites_);
       
        for(int j = 1; j < N_; ++j) {
            ampo += J*4.0,"Sz",j,"Sz",j+1;
            ampo += J*2.0,"S+",j,"S-",j+1;
            ampo += J*2.0,"S-",j,"S+",j+1;
        }
    
        if (BC == "PBC") {
            ampo += J*4.0,"Sz",1,"Sz",N_;
            ampo += J*2.0,"S+",1,"S-",N_;
            ampo += J*2.0,"S-",1,"S+",N_;
        }
    
        Hamiltonian_ = toMPO(ampo);
    }
    
    void XY1D(int J, string BC) {

        sites_ = SpinHalf(N_);
        auto ampo = AutoMPO(sites_);
    
        for(int j = 1; j < N_; ++j) {
            ampo += J*2.0,"S+",j,"S-",j+1;
            ampo += J*2.0,"S-",j,"S+",j+1;
        }
    
        if (BC == "PBC") {
            ampo += J*2.0,"S+",1,"S-",N_;
            ampo += J*2.0,"S-",1,"S+",N_;
        }
    
        Hamiltonian_ = toMPO(ampo);
    }
    
    void XXZ1D(double J, double K, string BC) {
        
        sites_ = SpinHalf(N_);
        auto ampo = AutoMPO(sites_);
    
        for(int j = 1; j < N_; ++j) {
            ampo += J*4.0,"Sz",j,"Sz",j+1;
            ampo += K*2.0,"S+",j,"S-",j+1;
            ampo += K*2.0,"S-",j,"S+",j+1;
        }
    
        if (BC == "PBC") {
            ampo += J*4.0,"Sz",1,"Sz",N_;
            ampo += K*2.0,"S+",1,"S-",N_;
            ampo += K*2.0,"S-",1,"S+",N_;
        }
    
        Hamiltonian_ = toMPO(ampo);
    }

    void InitializeState(string state_type) {
        // state_type options:
        // Neel, all up / down, random
    
        auto state = InitState(sites_);

        if (state_type == "neel") {
            uniform_int_distribution<int> distribution(0,1);
            
            for (int i = 1; i <= N_; ++i) {
                if (i%2 == 1) 
                    state.set(i,"Up");
                else
                    state.set(i, "Dn");
            }
        }

        if (state_type == "random") {
            uniform_int_distribution<int> distribution(0,1);
            mt19937 rd;
           
            for (int i = 1; i <= N_; ++i) {
                auto rng = distribution(rd);
                
                if (rng == 1)
                    state.set(i,"Up");
                else
                    state.set(i,"Dn");
    
            }
        }

        if (state_type == "Up") {
            for (int i = 1; i <= N_; ++i) {
                state.set(i,"Up");
            }
        }
        
        if (state_type == "Down") {
            for (int i = 1; i <= N_; ++i) {
                state.set(i,"Dn");
            }
        }

        psi0_ = MPS(state);
    }  

    ITensor ContractMPS(MPS psi) {
        ITensor psi_contract = psi.A(1);
        for (int i = 2; i <= N_; i++) {
            psi_contract *= psi.A(i);
        }
        return psi_contract;
    }
    
    void Run() {
        // TODO: maybe have the DMRG parameters user-input... not sure
        auto sweeps = Sweeps(100);
        sweeps.maxdim() = 50,50,100,100,200;
        
        auto [energy,psi] = dmrg(Hamiltonian_, psi0_, sweeps, "Quiet");
        psi_ = psi;
        printfln("\nGround State Energy = %.10f",energy/float(N_));
    }
};

#endif
