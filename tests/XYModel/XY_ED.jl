using LinearAlgebra
using DelimitedFiles

N = 2
J = -1

Dim = 2^N
Hamiltonian = zeros(Dim,Dim)

for Ket = 0:Dim-1
    for SpinIndex = 0:N-2
        bit1 = 2^SpinIndex
        bit2 = 2^(SpinIndex+1)
        
        # XX term
        Bra = (Ket ⊻ bit1) ⊻ bit2
        Hamiltonian[Bra+1,Ket+1] += 1

        ## YY term
        # check parity of spin chain
        parity = []
        for spin = 0:N-1
            # spin in binary (0,1)
            SpinVal = (Ket>>spin)&1
            push!(parity, SpinVal)
        end

        parity = sum(parity) % 2
        Hamiltonian[Bra+1,Ket+1] += (-1)^(parity+1)

    end
end

Diag = eigen(J*Hamiltonian)
GroundState = hcat(Diag.vectors[:, 1], zeros(Dim));
println(GroundState)
E0 = Diag.values[1] / N
E1 = Diag.values[2] / N

open("true_psi_N=$N", "w") do io
    writedlm(io, GroundState)
end

open("true_energies_N=$N", "w") do io
    writedlm(io, [E0, E1])
end
