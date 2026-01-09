module Trap

using IonSim
using QuantumOptics

export create_standard_chamber


    function create_standard_chamber(; B_field::Float64=4e-4)
        # We'll use Calcium-40, with the qubit encoded in the S1/2 and D5/2 levels.
        C = Ca40([("S1/2", -1/2, "g"),("D5/2", -1/2, "e")])

        chain = LinearChain(
            ions=[C],
            comfrequencies=(x=3e6, y=3e6, z=1e6),
            selectedmodes=(;z=[1]) # We only consider the axial mode of motion
        )

        # Define the laser that will drive the transition
        L = Laser()

        T = Chamber(
            iontrap=chain,
            B=4e-4,              # Magnetic field in Tesla
            Bhat=ẑ,              # Magnetic field direction
            lasers=[L]
        )
        # request are set here as properties of the laser object.
        # 'w' (frequency) is set by Δ (detuning) and λ (wavelength).
        # 'I' (intensity) is set by E (E-field amplitude).
        # 'phi' is set by L.ϕ.
        # 't' (time-dependence) is handled by the solver over a time span `tspan`.
        L.k = (x̂ + ẑ)/√2       # Laser wavevector
        L.ϵ = (x̂ - ẑ)/√2       # Laser polarization
        L.λ = transitionwavelength(C, ("g", "e"), T) # Set laser wavelength to be resonant with the transition
        L.Δ = 0.0
        pi_time = 2.5e-6
        E_amplitude = 1000000 #Efield_from_pi_time(pi_time, T, 1, 1, ("g", "e"))
        intensity!(L, t->E_amplitude) # This sets a constant (square) pulse

        return T
    end


end