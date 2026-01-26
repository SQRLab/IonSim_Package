module Pulse

using IonSim
using QuantumOptics

export simple_pulse, bfield_detune


    function simple_pulse(chamber; Duration::Float64=16)
        
        return h = experiment(chamber, Duration)
    end

    function pulse(T:: Chamber , tspan, pitime)
    # Define the laser that will drive the transition
    L = T.lasers[1]

    # Combine all components into a single Trap object, which represents the full experiment
    # This is the main object that holds the entire state of our physical system.
    
    pi2_time = pitime*1e6/2

    res_intensity = intensity_from_pitime(L, pitime, T.iontrap.ions[1], ("g", "e"), T)

    function intensity_funtion(t)
    if(t<=pi2_time)
        return res_intensity
    elseif(t>=tspan[end] - pi2_time)
        return res_intensity
    else
        return 0.0
    end

    end
    intensity!(L, intensity_funtion)

    function phase_funtion(t)
        if(t<=pi2_time)
            return 2*pi
        elseif(t>=tspan[end] - pi2_time)
            return pi
        else
            return 0.0
        end
    end

        phase!(L, phase_funtion)
        h = hamiltonian(T, timescale=1e-6, rwa_cutoff=Inf);
        return h
    end


    function experiment(T::Chamber, wait_time)
        pitime = 4e-6
        tspan = 0: 0.1: wait_time+4

        h = pulse(T, tspan, pitime)
        return h
    end

    function bfield_detune(T::Chamber, dB)
        bfield_fluctuation!(T, dB)
        return T
    end

end

