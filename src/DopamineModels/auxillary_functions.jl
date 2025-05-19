# external stimulus current (spike train or arbitrary function)
function I_stim(t; stim_start = 100, stim_end = 400, stim_amplitude = 200, pulse_width = 5)
    # 5 ms square pulses at 20 Hz between 100–400 ms
    if stim_start <= t <= stim_end && mod(t, 50) < pulse_width            # 20 Hz
        return stim_amplitude  # pA
    else
        return 0.0
    end
end