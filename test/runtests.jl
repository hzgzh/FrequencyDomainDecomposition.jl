using FrequencyDomainDecomposition
using Test

@testset "FrequencyDomainDecomposition.jl" begin
    # Write your tests here.
end

using MAT

y = matread("src/measurement1.mat")
y=y["y"]
fs =4096;                                   # Sample frequency
nfft=4096;                                  # Number of fft points
delta_freq = fs/nfft;                       # Frequency resolution
w = hanning(nfft);                          # Create window
sigs = MultiChannelSignals(y,fs)