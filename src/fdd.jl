



function decomposition(cpsd::CrossPowerSpectram,pos,n)
    data = cpsd.psd[:,:,pos]
    u,_,_ = svd(data)
    u[:,n]
end
"""
    mac(phi1,phi2)
    two mode similiar
"""
function mac(phi1,phi2)
    # This function calculates mac between phi1 and phi2
    mac= (abs(phi1'*phi2))^2/((phi1'*phi1)*(phi2'*phi2));
end

"""
    cpsd(sig::MultiChannelSignals,nfft;windows = hanning(nfft)) 
 Efficient cross spectral density estimation
 - `sig` : signal type of MultiChannelSignals
 - `nfft`: fft n points
 - `windows` : fft windows 
"""
function cpsd(sig::MultiChannelSignals,nfft;windows = hanning(nfft)) 
    #This algorithm uses the welch method with convolution and symmetry to
    #improve computational time of sd matrices compared to cpsd.m (without
    #spending time on other uneeded calculations)
    #*note 50% window overlap is implemented

    # define paramters
    
    fs = sample_rate(sig)
    NFFT=nfft ÷ 2 + 1
    f_vec=collect(LinRange(0,fs/2,NFFT))
    delta_freq = fs/nfft

    # Window parameters (50% overlap)
    y = sig.data
    nchannels = sig.nchannels
    len = sig.length

    win_size = nfft
    overlap = nfft ÷ 2
    step_size = win_size - overlap
    offset = 1:step_size:len - win_size+1
    m = length(offset)
       
    T = eltype(windows)
    # tic
    # preallocate space for the fft and psd results
    AP3  = zeros(Complex{T}, nchannels, nchannels, nfft)
    temp = zeros(Complex{T}, nfft, nchannels)
    ap = zeros(Complex{T},nfft)
    ap_conj = zeros(Complex{T},nfft)
    # loop through the data, one window at a time
    for i = 1:m
       for k = 1:nchannels
           # apply window
           temp[:,k] = y[offset[i]:offset[i]+win_size-1,k].*windows
           # calculate fft
           temp[:,k] = fft(temp[:,k])
       end
       # build lower triangular crosspower and sum blocks
       for nn=1:nchannels
           for mm=1:nn
               # convolution in frequency domain
               # producing a periodogram of "absolute value squared"
               ap .= temp[:,nn].*conj(temp[:,mm])
               AP3[nn,mm,:]=AP3[nn,mm,:]+ap[1:nfft]
                   # use hermetian symmetry of matrix to populate the
                   # upper triangular with the complex conjugate
                   if mm != nn
                       ap_conj .= conj(ap)
                       AP3[mm,nn,:]=AP3[mm,nn,:]+ap_conj[1:nfft]
                   end
           end
       end
    end

    # psd scaling factor (Brandt pg. 212)
    Sp = 1 / (nfft*delta_freq*sum(windows.^2))
    # scale and average the periodograms
    ap3 = (Sp/m)*AP3
    # select half spectra
    ap3 = ap3[:,:,1:nfft ÷ 2 + 1]
    # scale by 2 for redundatn nyquist frequencies except ignore DC and Nyquist value
    ap3[:,:,2:end-1] = ap3[:,:,2:end-1] * 2
    # t3 = toc;

    # plotting
    # for nn=1:2
    #     for mm=1:2
    #         h=reshape(ap3(nn,mm,:),1,NFFT);
    #         subplot(2,2,(nn-1)*2+mm)
    #         hold on; plot(f_vec,log10(abs(h)),'b');hold off;
    #     end
    # end 
    return CrossPowerSpectram(f_vec,ap3)
end

"""
    fdd(sig::MultiChannelSignals,w::Array{T,1}) where T
    Frequency Domain decomposition

** Argument **
- `y`: signal row observation column node
- `nfft`: fft points 
- `numbers_eigs`: selected numbers eigs at single frequecy point
- `windows`: fft windows hanning,hamming,cheb
"""
function fdd(sigs::MultiChannelSignals,nfft,numbers_eigs::Int64;windows = hanning(nfft))

    # *********************************************
    #               Measurement parameters 
    # *********************************************
    #fs = sig.sample_rate
    #nchannels = sig.nchannels
    #length = sig.length
    #y = sig.data
    #nfft = length(w)
    #delta_freq = fs/nfft   # Frequency resolution
    #NFFT=nfft ÷ 2 + 1
    

    # **************************************************************
    #               Power spectral density matrix 
    # **************************************************************
    #Gxy = zeros(sy[2],sy[2],NFFT)
    #F = zeros(sy[2],sy[2],NFFT)
    # efficient cross spectral density estimation
    p = cpsd(sigs,nfft;windows)
    Gxy = p.psd
    fvec = p.freq
    # **************************************************************
    #               Singular value decomposition 
    # **************************************************************
    
    sv = zeros(size(Gxy,3),numbers_eigs)
    for k = 1:size(Gxy,3)
        s = svdvals(Gxy[:,:,k])
            for i in 1:numbers_eigs
                sv[k,i] = s[i] # Singular values
            end
    end
    FDD(sigs,p,fvec,sv)
end

@recipe function f(sig::MultiChannelSignals,channel::Int64)
    @assert channel<nchannels(sig) "channel should not greater than $(nchannels(sig))"
    dt = 1/sig.sample_rate
    t = 0:dt:dt*(sig.length-1)
    legend := false
    xguide := "Time s"
    yguide := "channel $channel data"
    t,sig.data[:,channel]
end

@recipe function f(fdd::FDD)
    legend := false
    xguide := "frequencies Hz"
    yguide := "Amplitude Db"
    for i in 1:size(fdd.eigs)[2]
        @series begin
            fdd.freq,20log.(fdd.eigs[:,i])
        end
    end
end

@recipe function f(fdd::FDD,fp::Array{T,1},n::Int64) where T
    Δf = fdd.freq[end]/length(fdd.freq)
    nn = round.(fp./Δf,RoundUp)
    nn = Int.(nn)
    e  = fdd.eigs[nn]
    xguide := "frequencies Hz"
    yguide := "eigs db"
    legend := false
    @series begin
        fdd.freq,20log.(fdd.eigs[:,n])
    end
    
    @series begin
        seriestype := :scatter
        markershape := :circle
        marksershize := 4
        fp,20log.(e)
    end
end

@recipe function f(fdd::FDD,freq::Float64,n::Int64)
    xguide := "nodes"
    yguide := "Amplitude"
    abs.(peak_modal(fdd,freq,n))
end
"""
    findpeak(fdd::FDD,n::Int64;σ=2,threshold=10)
    find spectram peak position 
- `n`: eigs numbers
- `σ`: peak width
- `threshold`: peak should greater than threshold 0-100(from min to max)
"""
function findpeak(fdd::FDD,n::Int64;σ=2,threshold=10)
    data = fdd.eigs

    @assert n <= size(data)[2] "n should not greater than $(size(data)[2])"
    _,fp = peakfinder(20log.(data[:,n]);σ,threshold)
    return fp
end


function findmaxeigs(fdd::FDD,n::Int64;σ=2,threshold=10) 
    fp = findpeak(fdd,n;σ,threshold)
    Δf = fdd.freq[end]/length(fdd.freq)
    nn = round.(fp./Δf,RoundUp)
    nn = Int.(nn)
    fp[sortperm(fdd.eigs[nn,n],rev=true)],sort(fdd.eigs[nn,n],rev=true)
end

"""
    peak_modal(fdd::FDD,fp::Array{T,1},n::Int64) where T
    node modal
** Arguments **
- `fdd`: type of FDD
- `fp` : peak position
"""
function peak_modal(fdd::FDD,fp::Array{T,1},n::Int64) where T
    Δf = fdd.freq[end]/length(fdd.freq)
    nn = round.(fp./Δf,RoundUp)
    nn = Int.(nn)
    ns = nchannels(fdd.signals)
    v = zeros(Complex{T},ns,length(nn))
    for i in 1:length(nn)
        v[:,i] = decomposition(fdd.psd,nn[i],n)
    end
    fp,v
end

function peak_modal(fdd::FDD,fp::Float64,n::Int64)
    Δf = fdd.freq[end]/length(fdd.freq)
    nn = round.(fp/Δf,RoundUp)
    nn = Int(nn)
    decomposition(fdd.psd,nn,n)
end

export cpsd,fdd,peak_modal,findmaxeigs,f,mac


