module Simulation

using IonSim
using QuantumOptics

export run_simulation, run_simulation_hamiltonian

    function run_simulation(chamber, tspan)
        
        
        mode = zmodes(chamber)[1]
        #mode.N = 10
        C = Ca40([("S1/2", -1/2, "g"),("D5/2", -1/2, "e")])

        # Construct the time-independent Hamiltonian for the system
        # For a simple square pulse, the Hamiltonian is not time-dependent.
        # We set a high rwa_cutoff to ensure all terms are included.
        H = hamiltonian(chamber, rwa_cutoff=Inf)
        ψ_mode = fockstate(mode[1].basis, 0)
        ψ₀ = C["g"] ⊗ ψ_mode

        # Solve the time evolution using the Schrödinger equation solver from QuantumOptics.jl
        tout, sol = timeevolution.schroedinger_dynamic(tspan, ψ₀, H)

        # 5. Analyze and Visualize the Results
        # Calculate the population in the excited state |e⟩ over time
        excited_pop = expect(ionprojector(chamber, "e"), sol)
        return tout, excited_pop

    end

    function run_simulation_hamiltonian(t, h, tspan)
         mode = zmodes(t)[1]
        C = Ca40([("S1/2", -1/2, "g"),("D5/2", -1/2, "e")])
        H = hamiltonian(t, rwa_cutoff=Inf)
        ψ_mode = fockstate(mode[1].basis, 0)
        ψ₀ = C["g"] ⊗ ψ_mode
        tout, sol = timeevolution.schroedinger_dynamic(tspan, ψ₀, h)
        # ψ_mode = fockstate(mode[1].basis, 0)
        # ψ₀ = C["g"] ⊗ ψ_mode
        excited_pop = expect(ionprojector(t, "e"), sol)
        return tout, excited_pop
    end

end