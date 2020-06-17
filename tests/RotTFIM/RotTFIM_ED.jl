using LinearAlgebra
using DelimitedFiles

N = 2
J = 0.5
h = 3.0

Dim = Int(2^N)
Hamiltonian = zeros(Dim,Dim)

for Ket = 0:Dim-1
    diagonal = 0.
    for SpinIndex = 0:N-1
        Spin = 2*((Ket>>SpinIndex)&1) - 1
        diagonal += Spin
    end
    # fills entries the other way for some reason...
    Hamiltonian[Dim - (Ket), Dim - (Ket)] = -h*diagonal

    for SpinIndex = 0:N-2
        Spin1 = 2*((Ket>>SpinIndex)&1) - 1
        NextIndex = SpinIndex + 1
        Spin2 = 2*((Ket>>NextIndex)&1) - 1

        bit1 = 2^SpinIndex
        bit2 = 2^NextIndex
        Bra = Ket ⊻ bit1 ⊻ bit2
        Hamiltonian[Bra+1,Ket+1] = -J
    end
end

Diag = eigen(Hamiltonian)

GroundState = hcat(Diag.vectors[:, 1], zeros(Dim));
E0 = Diag.values[1] / N
E1 = Diag.values[2] / N

open("true_psi_N=$N"*"_J=$J"*"_h=$h", "w") do io
    writedlm(io, GroundState)
end

open("true_energies_N=$N"*"_J=$J"*"_h=$h", "w") do io
    writedlm(io, [E0, E1])
end
