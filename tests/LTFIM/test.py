from qucumber.utils import unitaries
import numpy as np
import torch
import matplotlib.pyplot as plt
import sys

sys.path.insert(1, '../')

import test_utils

'''
- Generate samples from DMRG, ED wavefunctions
- Bin these samples along with the samples from the MPS algorithm

Do this for all bases
'''

unitary_dict = {k: test_utils.convert_torch_cplx(
    v) for k, v in unitaries.create_dict().items()}

N = 2

bases = np.array([
    ["Z", "Z"],
    ["X", "X"],
    ["Y", "Y"],
    ["Z", "X"],
    ["X", "Z"],
    ["Z", "Y"],
    ["Y", "Z"],
    ["Y", "X"],
    ["X", "Y"]
])

ED_psi = np.loadtxt("true_psi_N=2_J=0.5_h=1.0_Omega=3.0")
ED_psi = ED_psi[:, 0] + 1j*ED_psi[:, 1]

num_samples = 10000
DMRG_samples = np.loadtxt("MPS_sampling_algo_N=2_J=0.5_h=1_Omega=3", dtype=int)

for i in range(len(bases)):
    basis = bases[i]
    print(basis)
    psi_r = test_utils.rotate_psi(unitary_dict, basis, ED_psi)
    probs = (psi_r * np.conj(psi_r)).real
    print(psi_r)
    probs = probs[probs > 1e-13]
    print(probs)

    ED_inds, ED_samples = test_utils.gen_samples(num_samples, N, psi_r)

    ALGO_samples = DMRG_samples[10000*i:10000*(i+1)]
    ALGO_inds = test_utils.gen_inds_from_samples(ALGO_samples)

    ED_uniques, ED_counts = np.unique(ED_inds, return_counts=True)
    ALGO_uniques, ALGO_counts = np.unique(ALGO_inds, return_counts=True)
    ED_counts = ED_counts / len(ED_inds)
    ALGO_counts = ALGO_counts / len(ALGO_inds)

    plt.figure()
    plt.bar(ED_uniques+0.25, ED_counts, color='blue',
            label="ED samples", align='center', width=0.25)
    plt.bar(ED_uniques, probs, color='black',
            label="True prob", align='center', width=0.25)
    plt.bar(ALGO_uniques-0.25, ALGO_counts, color='green',
            label="ALGO samples", align='center', width=0.25)
    plt.title(basis[0]+basis[1])
    plt.xlabel("Basis state index")
    plt.ylabel("Fractional frequency")
    plt.legend()
    plt.xticks([0, 1, 2, 3])
    plt.savefig("Basis={}.pdf".format(
        basis[0]+basis[1]), dpi=500, bbox_inches='tight')
