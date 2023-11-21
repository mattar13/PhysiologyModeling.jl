#Auxillary equations not needing parameters
αmKV(v) = (5 * (100 - v)) / (exp((100 - v) / 42) - 1)
βmKV(v) = 9 * exp(-(v - 20) / 40)
αhKV(v) = 0.15 * exp(-v / 22)
βhKV(v) = 0.4125 / (exp((10 - v) / 7) + 1)