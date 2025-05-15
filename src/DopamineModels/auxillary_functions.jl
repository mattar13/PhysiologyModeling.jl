# external stimulus current (spike train or arbitrary function)
function I_stim(t)
    # 5 ms square pulses at 20 Hz between 100–400 ms
    if 100 <= t <= 400 && mod(t, 50) < 5            # 20 Hz
        return 200.0  # pA
    else
        return 0.0
    end
end