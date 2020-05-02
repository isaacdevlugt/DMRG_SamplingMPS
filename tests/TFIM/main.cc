#include "itensor/all.h"
#include <cstdlib>
#include <sstream>
#include <random>
#include <fstream>
#include <iostream>

#include <string>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../../dmrg.h"
#include "../../sampler.h"

using namespace std;
using namespace itensor;

int main() {

    int N = 2;
    int num_samples = 10000;
    string bc = "OBC";
    ifstream bases_file("../../bases.txt");
    
    DMRG dmrg_(N);
    // TFIM at critical point
    double J = 0.5;
    double h = 3.0;
    dmrg_.TFIM1D(J, h, bc);
    dmrg_.InitializeState("random");
    dmrg_.Run();

    MPS psi;
    psi = dmrg_.GetPsi();
    SiteSet sites = dmrg_.GetSiteSet();
  
    // save the DMRG wavefunction 
    ITensor DMRG_psi = dmrg_.ContractMPS(psi);
    auto Inds = DMRG_psi.inds();
    stringstream DMRG_psi_path;
    DMRG_psi_path << "DMRG_psi_N=" << N << "_J=" << J << "_h=" << h << ends;
    ofstream DMRG_psi_file(DMRG_psi_path.str());

    for (int n1 = 1; n1 <= dim(Inds[0]); n1++) {
        for (int n2 = 1; n2 <= dim(Inds[1]); n2++) {
            DMRG_psi_file << DMRG_psi.cplx(Inds[0](n1), Inds[1](n2)) << endl;
        }
    }

    stringstream sample_path_;
    stringstream sample_bases_path_;
    sample_path_ << "MPS_sampling_algo_N=" << N << "_J=" << J << "_h=" << h << ends;
    sample_bases_path_ << "MPS_sampling_algo_bases_N=" << N << "_J=" << J << "_h=" << h << ends;
    auto sample_path = sample_path_.str();
    auto sample_bases_path = sample_bases_path_.str();

    Sampler sampler(N, psi, sites, num_samples, sample_path, sample_bases_path);
    sampler.LoadBases(bases_file);
    sampler.Sample();
}
