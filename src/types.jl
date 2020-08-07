import Base:length,show
"""
    MultiChannelSignals{T<:AbstractFloat,G<:Real}
    multichannels signal row represent measurements,column represent node

- `data`       :  signals data
- `nchannels`  :  channels
- `length`     :  signal length
- `sample_rate`:  signal sample rate
"""
struct MultiChannelSignals{T<:AbstractFloat,G<:Real}
    data::Array{T,2}
    nchannels::Int64
    length::Int64
    sample_rate::G
    
    function MultiChannelSignals(sigs::Array{T,2},fs::G) where {T,G}
        len,nchannels = size(sigs)
        @assert len>nchannels "row should be observation,column should be nodes,rows should great than cols"
        @assert nchannels>2 "signal channels should great than 1"
        return new{T,G}(sigs,nchannels,len,fs)
    end
end

function Base.show(io::IO,signal::MultiChannelSignals)
    println("signal channels : $(nchannels(signal))")
    println("signal length   : $(length(signal))")
    println("signal sample rate : $(sample_rate(signal))")
end

"""
    nchannels(sig::MultiChannelSignals)
    signal channels
"""
nchannels(sig::MultiChannelSignals)=sig.nchannels
"""
    Base.length(sig::MultiChannelSignals)
    every channel signal length
"""
Base.length(sig::MultiChannelSignals)=sig.length
"""
    sample_rate(sig::MultiChannelSignals)
    every channel signal sample rate
"""
sample_rate(sig::MultiChannelSignals)=sig.sample_rate
"""
    channel(sig::MultiChannelSignals,channel::Int64)
    single channel signal data
"""
channel(sig::MultiChannelSignals,channel::Int64) = sig.data[:,channel]

"""
    CrossPowerSpectram{G,T}
    cross power spectram of multichannel signal
- `freq`: frequency of spectram
- `psd` : cross power spectgram 
"""
struct CrossPowerSpectram{G,T}
    freq::Array{G,1}
    psd::Array{T,3}
    function CrossPowerSpectram(freq::Array{G,1},psd::Array{T,3}) where {G,T}
        new{G,T}(freq,psd)
    end
end


"""
    FDD{T}
    result of frequency domain decomposition
- `signal`: type of MultiChannelSignals
- `psd`   : type of CrossPowerSpectram
- `freq`  : frequencies vector
- `eigs`  : svd eigs
"""
struct FDD{T,G}
    signals::MultiChannelSignals
    psd::CrossPowerSpectram
    freq::Array{T,1}
    eigs::Array{G,2}
    function FDD(signals,psd,freq::Array{T,1},eigs::Array{G,2}) where {T,G}
        new{T,G}(signals,psd,freq,eigs)
    end
end



export MultiChannelSignals,FDD,CrossPowerSpectram
export nchannels,sample_rate,channel