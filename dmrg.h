#ifndef DMRG_H
#define DMRG_H

#include "itensor/all.h"
#include <cstdlib>
#include <random>
#include <fstream>

using namespace std;
using namespace itensor;

class DMRG
{

    int N_;
    MPO Hamiltonian_;
    SiteSet sites_;
    MPS psi_;
    double energy_;
    double abs_magnetization_;
    MPS psi0_;

public:
    DMRG(int N) : N_(N) {}

    inline MPS GetPsi()
    {
        return psi_;
    }

    inline double GetEnergy()
    {
        return energy_;
    }

    inline double GetAbsM()
    {
        return abs_magnetization_;
    }

    inline SiteSet GetSiteSet()
    {
        return sites_;
    }

    void Rydberg(double R, double h, double Omega)
    {
        vector<vector<double>> Vij;

        // inefficient way to store Vij since a lot of repitition
        for (int i = 1; i < N_; ++i)
        {
            vector<double> tmp;
            
            for (int j = i + 1; j <= N_; ++j)
            {
                auto vij = 1/pow(R*(j-i), 6.0);
                tmp.push_back(vij);
            }
            
            Vij.push_back(tmp); 
        }

        sites_ = SpinHalf(N_, {"ConserveQNs=", false});
        auto ampo = AutoMPO(sites_);

        for (int i = 1; i < N_; ++i)
        {
            for (int j = i + 1; j <= N_; ++j)
            {
                ampo += -Vij[i-1][j-i-1] * 4.0, "Sz", i, "Sz", j;
                ampo += -Vij[i-1][j-i-1] * 2.0, "Sz", i;
                ampo += -Vij[i-1][j-i-1] * 2.0, "Sz", j;
            }

            ampo += -h * 2.0, "Sz", i;
            ampo += -Omega * 2.0, "Sx", i;
        }

        ampo += -h * 2.0, "Sz", N_;
        ampo += -Omega * 2.0, "Sx", N_;

        Hamiltonian_ = toMPO(ampo); 
    }

    void LTFIM1D(double J, double h, double Omega, string BC)
    {
        sites_ = SpinHalf(N_, {"ConserveQNs=", false});
        auto ampo = AutoMPO(sites_);

        for (int j = 1; j < N_; ++j)
        {
            ampo += -J * 4.0, "Sz", j, "Sz", j + 1;
            ampo += -h * 2.0, "Sx", j;
            ampo += -Omega * 2.0, "Sz", j;
        }
        ampo += -h * 2.0, "Sx", N_;
        ampo += -Omega * 2.0, "Sz", N_;

        if (BC == "PBC")
        {
            ampo += -J * 4.0, "Sz", 1, "Sz", N_;
        }

        Hamiltonian_ = toMPO(ampo);
    }

    void RotTFIM1D(double h, string BC)
    {

        sites_ = SpinHalf(N_, {"ConserveQNs=", false});
        auto ampo = AutoMPO(sites_);

        for (int j = 1; j < N_; ++j)
        {
            ampo += -4.0, "Sx", j, "Sx", j + 1;
            ampo += -h * 2.0, "Sz", j;
        }
        ampo += -h * 2.0, "Sz", N_;

        if (BC == "PBC")
        {
            ampo += -4.0, "Sx", 1, "Sx", N_;
        }

        Hamiltonian_ = toMPO(ampo);
    }

    void TFIM1D(double h, string BC)
    {

        sites_ = SpinHalf(N_, {"ConserveQNs=", false});
        auto ampo = AutoMPO(sites_);

        for (int j = 1; j < N_; ++j)
        {
            ampo += -4.0, "Sz", j, "Sz", j + 1;
            ampo += -h * 2.0, "Sx", j;
        }
        ampo += -h * 2.0, "Sx", N_;

        if (BC == "PBC")
        {
            ampo += -4.0, "Sz", 1, "Sz", N_;
        }

        Hamiltonian_ = toMPO(ampo);
    }

    void Heisenberg1D(double J, string BC)
    {

        sites_ = SpinHalf(N_, {"ConserveQNs=", false});
        auto ampo = AutoMPO(sites_);

        for (int j = 1; j < N_; ++j)
        {
            ampo += J * 4.0, "Sz", j, "Sz", j + 1;
            ampo += J * 2.0, "S+", j, "S-", j + 1;
            ampo += J * 2.0, "S-", j, "S+", j + 1;
        }

        if (BC == "PBC")
        {
            ampo += J * 4.0, "Sz", 1, "Sz", N_;
            ampo += J * 2.0, "S+", 1, "S-", N_;
            ampo += J * 2.0, "S-", 1, "S+", N_;
        }

        Hamiltonian_ = toMPO(ampo);
    }

    void XY1D(int J, string BC)
    {

        //auto q = QN({"Sz",0})
        // TODO: should have Sz=0 conserved here but there is a bug in sampler.h
        // --> sites_ = SpinHalf(N_, {"ConserveQNs=", true});
        // bug only happens when wanting to sample in other bases
        // leave true if only Z basis samples

        sites_ = SpinHalf(N_, {"ConserveQNs=", true});
        auto ampo = AutoMPO(sites_);

        for (int j = 1; j < N_; ++j)
        {
            ampo += J * 2.0, "S+", j, "S-", j + 1;
            ampo += J * 2.0, "S-", j, "S+", j + 1;
        }

        if (BC == "PBC")
        {
            ampo += J * 2.0, "S+", 1, "S-", N_;
            ampo += J * 2.0, "S-", 1, "S+", N_;
        }

        Hamiltonian_ = toMPO(ampo);
    }

    void XXZ1D(double J, double K, string BC)
    {

        sites_ = SpinHalf(N_);
        auto ampo = AutoMPO(sites_);

        for (int j = 1; j < N_; ++j)
        {
            ampo += J * 4.0, "Sz", j, "Sz", j + 1;
            ampo += K * 2.0, "S+", j, "S-", j + 1;
            ampo += K * 2.0, "S-", j, "S+", j + 1;
        }

        if (BC == "PBC")
        {
            ampo += J * 4.0, "Sz", 1, "Sz", N_;
            ampo += K * 2.0, "S+", 1, "S-", N_;
            ampo += K * 2.0, "S-", 1, "S+", N_;
        }

        Hamiltonian_ = toMPO(ampo);
    }

    void InitializeState(string state_type)
    {
        // state_type options:
        // Neel, all up / down, random

        auto state = InitState(sites_);

        if (state_type == "neel_Sz=0")
        {
            for (int i = 1; i <= N_; ++i)
            {
                if (i % 2 == 1)
                    state.set(i, "Up");
                else
                    state.set(i, "Dn");
            }
                
        }

        if (state_type == "neel")
        {
            for (int i = 1; i <= N_; ++i)
            {
                if (i % 2 == 1)
                    state.set(i, "Up");
                else
                    state.set(i, "Dn");
            }
        }

        if (state_type == "random")
        {
            uniform_int_distribution<int> distribution(0, 1);
            mt19937 rd;

            for (int i = 1; i <= N_; ++i)
            {
                auto rng = distribution(rd);

                if (rng == 1)
                    state.set(i, "Up");
                else
                    state.set(i, "Dn");
            }
        }

        if (state_type == "Up")
        {
            for (int i = 1; i <= N_; ++i)
            {
                state.set(i, "Up");
            }
        }

        if (state_type == "Down")
        {
            for (int i = 1; i <= N_; ++i)
            {
                state.set(i, "Dn");
            }
        }

        psi0_ = MPS(state);
    }

    ITensor ContractMPS(MPS psi)
    {
        ITensor psi_contract = psi.A(1);
        for (int i = 2; i <= N_; i++)
        {
            psi_contract *= psi.A(i);
        }
        return psi_contract;
    }

    void Run()
    {
        auto sweeps = Sweeps(100);
        sweeps.maxdim() = 50, 50, 100, 100, 200;

        auto [energy, psi] = dmrg(Hamiltonian_, psi0_, sweeps, "Quiet");
        psi_ = psi;
        //sites_ = GetSiteSet();
        
        // measure magnetization
        Real Sz;
        for (int j = 1; j <= N; ++j) 
        {
            psi_.position(j);
            Real Szj = elt(psi_(j)
                       * op(sites_, "Sz", j)
                       * dag(prime(psi_(j), "Site")));
            Sz += Szj;
        }

        abs_magnetization_ = abs(Sz) / float(N_);
        energy_ = energy / float(N_);
        printfln("\nGround State E = %.10f", energy_);
        printfln("\nGround State |M| = %.10f", abs_magnetization_);
    }
};

#endif
