"""
    wTSM(y::Array{T,1}, yfit::Array{T,1}, w::Array{T,1}; 
        iter::Integer = 2, nptperyear::Integer = 46, wfact::Float64 = 2.0) 

Weight updating method in TIMESAT

# Author
Translated from TIMESAT to Julia by Dongdong Kong (20200808)
"""
function wTSM(y::AbstractArray{T,1}, yfit::AbstractArray{T,1}, w::AbstractArray{T,1};
    iter::Integer = 2, nptperyear::Integer = 46, wfact::Float64 = 0.5) where {T <: AbstractFloat}

    n = length(y)
    m = sum(w .> 0.5)
    w_ceil = ceil.(w)
    wnew = deepcopy(w)

    # yfit = y
    yfitmean = sum(yfit .* w_ceil / m)
    yfitstd = sqrt(sum(@. ((yfit - yfitmean) * w_ceil )^2)/(m-1))
    deltaT = fld(nptperyear, 7)
    
    @inbounds for i = 1:n
        m1 = max(1, i - deltaT);
        m2 = min(n, i + deltaT);
        idx = m1:m2
        # println(idx)
        yi = yfit[idx]
        yi_min = minimum(yi)
        yi_max = maximum(yi)

        # Adjust the weights dependent on if the values are above or below the 
        # fitted values
        if (y[i] < yfit[i] - 1e-8)
            #  if (yi_min > yfitmean){
            if (yi_min > yfitmean || iter < 2)
                # If there is a low variation in an interval, i.e. if the interval
                # is at a peak or at a minima compute the normalized distance
                # between the data point and the fitted point.
                if (yi_max - yi_min < 0.8 * yfitstd) 
                    ydiff = 2 *(yfit[i] - y[i])/yfitstd;
                else
                    ydiff = 0;
                end
                # Use the computed distance to modify the weight. Large distance
                # will give a small weight
                wnew[i] = wfact * w[i] * exp(-ydiff*ydiff);
            end
        end
    end   
    wnew 
end
