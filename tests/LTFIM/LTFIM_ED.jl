using LinearAlgebra
using DelimitedFiles

N = 2
J = 0.5
h = 1.0
Ω = 3.0

Dim = 2^N
Hamiltonian = zeros(Dim,Dim)

for Ket = 0:Dim-1
    diagonal = 0.
    for SpinIndex = 0:N-2
        Spin1 = -(2*((Ket>>SpinIndex)&1) - 1)
        NextIndex = SpinIndex + 1
        Spin2 = -(2*((Ket>>NextIndex)&1) - 1)
        diagonal = diagonal - J*Spin1*Spin2
    end

    for SpinIndex = 0:N-1
        Spin = -(2*((Ket>>SpinIndex)&1) - 1)
        diagonal = diagonal - Ω*Spin
    end

    Hamiltonian[Ket+1,Ket+1] = diagonal

    for SpinIndex = 0:N-1
        bit = 2^SpinIndex
        Bra = Ket ⊻ bit
        Hamiltonian[Bra+1,Ket+1] = -h
    end
end

Diag = eigen(Hamiltonian)
GroundState = hcat(Diag.vectors[:, 1], zeros(Dim));
E0 = Diag.values[1] / N
E1 = Diag.values[2] / N

global magnetization = 0
for Ket = 0:Dim-1  #Loop over Hilbert Space
    SumSz = 0.
    for SpinIndex = 0:N-1  #Loop over spin index (base zero, stop one spin before the end of the chain)
        Spin1 = -(2*((Ket>>SpinIndex)&1) - 1)
        SumSz += Spin1 #spin is +1 or -1
    end
    global magnetization += SumSz*GroundState[Ket+1]^2  #Don't forgot to square the coefficients...
end

global magnetization = abs(magnetization) / N

open("true_psi_N=$N"*"_J=$J"*"_h=$h"*"_Omega=$Ω", "w") do io
    writedlm(io, GroundState)
end

open("true_energies_absmag_N=$N"*"_J=$J"*"_h=$h"*"_Omega=$Ω", "w") do io
    writedlm(io, [E0, E1, magnetization])
end
