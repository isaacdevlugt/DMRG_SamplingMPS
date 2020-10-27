using LinearAlgebra
using DelimitedFiles

N = 4
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

        # YY term
        # check parity of spin chain
        tmp = Ket & Bra
        parity = 0
        for spin = 0:N-1
            # spin in binary (0,1)
            SpinVal = (tmp>>spin)&1
            parity += SpinVal
        end
    
        #parity = sum(

        #parity = sum(parity) # really mod 2... but doesn't matter
        Hamiltonian[Bra+1,Ket+1] += (-1)^(parity+1)

    end
end

for i=1:2^N
    println(Hamiltonian[i,:])
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
