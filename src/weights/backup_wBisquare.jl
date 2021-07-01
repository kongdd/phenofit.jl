using Statistics

"""
wBisquare weights updating 

wBisquare(y, yfit, w; iter = 2, wmin, to_upper = true)

# arguments
- `iter`: not used

# Bad points
- 1. under the yfit, in the growing season (yfit > 0.3 * A + ymin)
- 2. 
"""
function wBisquare(y::Array{T,1}, yfit::Array{T,1}, w::Array{T2,1};
    # QC_flag;
    iter::Integer = 2, 
    wmin::Float64 = 0.05, 
    step::Float64 = 0.5, 
    to_upper = true) where {T <: AbstractFloat, T2 <: AbstractFloat}
    
    ymax = maximum(yfit)
    ymin = minimum(yfit)
    A = ymax - ymin;

    re     = yfit .- y
    re_abs = abs.(re)
    sc     = 6*median(re_abs)

    trs_high = 0.7
    trs_low  = 0.4

    # 最保险的方法，获取每年的ylu，然后判断是ingrowing or ungrowing
    # println("threshold：high=", trs_high * A + ymin, ", low=", trs_low * A + ymin)
    I_bad_high = @.( ((re > 0) & (yfit > trs_high * A + ymin)) ) # middle of GS, upper envelope
    I_bad_low = @.( ((re < 0) & (yfit < trs_low * A + ymin)) )
    I_bad = I_bad_high .| I_bad_low

    ## 1.1 坏的点，一窝端。要坏一起坏，重新洗牌，分配权重
    wnew = w; 
    # multiply 0.5 make sure bad values have low weight
    wnew[I_bad] = @.( (1 - (re_abs[I_bad]/sc)^2)^2 * 0.5) # - step
    # wnew2 = @.( (1 - (re_abs[I_bad]/sc)^2)^2 ); # template
    # wnew[I_bad] = wnew[I_bad] .* w[I_bad]

    ## 1.2 采取一些抢救措施，拯救坏的点
    # y[I_bad_high] .= yfit[I_bad_high]
    y[I_bad] .= yfit[I_bad]
    
    ## 2 好的点，予以奖励，增加权重
    I_good = @.( ((re < 0) & (yfit > trs_high * A + ymin)) | 
                 ((re > 0) & (yfit < trs_low * A + ymin)))
    wnew[I_good] = wnew[I_good] .+ step

    ## 3 异常值赋予最低权重，改写对应的value
    I_outlier = (re .> sc) # .& (QC_flag .!== Int(1)
    wnew[I_outlier] .= wmin
    wnew[wnew .< wmin] .= wmin
    y[I_outlier] = yfit[I_outlier]
    
    ## also need to change y to improve the upper envelope performance
    # println("inside wBisquare: ", sum(y))
    # wmax = 2.0
    # wnew[wnew .>= 1.0] .= 1.0
    wnew
end
