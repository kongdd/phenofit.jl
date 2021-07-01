"""
    findpeaks(x::AbstractVector{T};
        nups::Int=1,
        ndowns::Int=nups,
        zerostr::Char='0',
        peakpat=nothing, 
        minpeakheight=typemin(T), 
        minpeakdistance::Int=1,
        A_max=zero(T),
        A_min=zero(T),
        npeaks::Int=0,
        sortstr=false)


# Filter condition
```julia
p.value >= minpeakheight && 
    p.value - max(x[p.start], x[p.stop]) >= A_min && 
    p.value - min(x[p.start], x[p.stop]) >= A_max 
```

# Examples
```
findpeaks(x; nups, ndowns, zerostr, peakpat=nothing, minpeakheight=typemin(T), threshold=zero(T), 
    npeaks=0, sortstr=false)
```

# References
@author: Gerhard Aigner   
https://github.com/halleysfifthinc/Peaks.jl/issues/11#issuecomment-689998279
"""
function findpeaks(x::AbstractVector{T};
    nups::Int=1,
    ndowns::Int=nups,
    zerostr::Char='0',
    peakpat=nothing, 
    options...) where {T <: Real}
    
    zerostr ∉ ('0', '+', '-') && error("zero must be one of `0`, `-` or `+`")

    # generate the peak pattern with no of ups and downs or use provided one
    peakpat = Regex(peakpat === nothing ? "[+]{$nups,}[-]{$ndowns,}" : peakpat)

    # transform x into a "+-+...-+-" character string
    xs = String(map(diff(x)) do e
        e < 0 && return '-'
        e > 0 && return '+'
        return zerostr
    end)

    grps = filter(x -> length(x) > 0, findall(peakpat, xs))
    # find index positions and maximum values
    peaks = map(grps) do m
        v, i = findmax(@view x[m]) 
        # fix extreme value positions on the plateau
        i = floor(Int, median(findall(@view(x[m]) .== v)))
        start = first(m)
        stop = last(m)+1
        idx = first(m) + i - 1
        (;start=start, idx=idx, stop=stop, 
            val_start = x[start], val=v, val_stop = x[stop], 
            diff_min = minimum([v - x[start], v - x[stop]]), 
            diff_max = maximum([v - x[start], v - x[stop]]))
    end

    if (length(options) > 0); 
        return filter_peaks(x, peaks; options...)
    else
        return peaks
    end
end

"""
    filter_peaks(peaks, 
        minpeakheight=typemin(T), minpeakdistance::Int=1,
        A_max=zero(T), A_min=zero(T), 
        npeaks::Int=0,
        sortstr=false, 
        history=false,
        ignored...)

# Parameters
- `history`: if true, removed reason will be returned

"""
function filter_peaks(x::AbstractVector{T}, peaks;
    minpeakheight=typemin(T), minpeakdistance::Int=1,
    A_max=0, A_min=0, 
    npeaks::Int=0,
    sortstr=false, 
    history=false,
    verbose=true,
    ignored...) where T <: Real

    if (verbose); println("A_max = $A_max, A_min = $A_min"); end
    # 记录每个点被剔除的原因
    # eliminate peaks that are too low
    peaks = deepcopy(peaks)
    n = length(peaks)

    d = peaks |> DataFrame
    removal = trues(n) # if true, then will be removed
    status = Array{String}(undef, n)
    for i = 1:length(peaks)
        p = peaks[i]
        if p.val < minpeakheight
            status[i] = "minpeakheight"
        elseif abs(p.diff_min) < A_min
            # println(abs(p.val - max(x[p.start], x[p.stop])))
            status[i] = "A_min"
        elseif abs(p.diff_max) < A_max
            status[i] = "A_max"
        else
            status[i] = ""
            removal[i] = false
        end
    end
    deleteat!(peaks, removal)
    # filter!(peaks) do p
    #     # p.value >= minpeakheight && p.value - max(x[p.start], x[p.stop]) >= threshold
    #     p.val >= minpeakheight && 
    #         p.val - max(x[p.start], x[p.stop]) >= A_min && 
    #         p.val - min(x[p.start], x[p.stop]) >= A_max 
    # end

    # sort according to peak height
    if sortstr || minpeakdistance > 1
        sort!(peaks, by=x -> x.val; rev=true)
    end
    
    # find peaks sufficiently distant
    all = d[!, :idx]
    left_1 = map(x -> x.idx, peaks)

    if minpeakdistance > 1
        removal = falses(length(peaks))
        for i in 1:length(peaks)
            removal[i] && continue
            for j in 1:length(peaks)
                removal[j] && continue
                dist = abs(peaks[i].idx - peaks[j].idx)
                removal[j] = 0 < dist < minpeakdistance 
            end
        end
        println(removal)
        deleteat!(peaks, removal)
    end
    left_2 = map(x -> x.idx, peaks)
    bads = setdiff(left_2, left_1)
    
    if (!is_empty(bads))
        ind_bad = deleteat!(collect(1:n), indexin(bads, all)) # look for deleted
        status[ind_bad] .= "minpeakdistance";        
    end

    # @show peaks left all status
    if history
        d[:, :status] = status
        d
    else
        # if (is_empty(left)); return nothing; end
        npeaks > 0 && resize!(peaks, min(length(peaks), npeaks))
        peaks
    end
    # Return only the first 'npeaks' peaks
    # history ? d : peaks # return
end

export filter_peaks
