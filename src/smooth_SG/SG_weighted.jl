using LinearAlgebra
# use class object here will improve SG performance significantly
Base.@kwdef mutable struct mod_SG{FT}
    # n    ::Integer
    # halfwin::Integer
    # FT = Float64
    # y :: AbstractArray{FT, 1}
    # w :: AbstractArray{FT2, 1}
    S    ::AbstractArray{Int, 2} # [frame, d+1], static array
    SMat ::AbstractArray{FT, 2}   # [frame, d+1], temp variable, update everytime
    T    ::AbstractArray{FT, 2}   # [d+1  , frame]
    B    ::AbstractArray{FT, 2}   # [frame, frame]
end

function init_mSG(halfwin = 1, d = 2)
    frame = halfwin*2 + 1;

    FT = Float64
    S    = zeros(Int, frame, d+1)
    SMat = zeros(FT, frame, d+1)
    T    = zeros(FT, frame, d+1)
    B    = zeros(FT, frame, frame)
    mod_SG(S, SMat, T, B)
end

function sgmat_wB!(w::AbstractArray{T1}, m::mod_SG) where T1 <: Real
    multiply_w_sqrt!(m.S, w, m.SMat)
    r = qr(m.SMat) # double
    m.T = as_SMatrix(r.R') \ m.S'; # four times faster
    
    mul!(m.B, m.T', m.T)
    # m.B = m.T' * m.T;
    multiply_col!(m.B, w)
end

"""
    wSG(y::Array{T, 1}, w::Array{T2, 1}; halfwin=1, d=2)

weighted Savitzky Golay filter
"""
function wSG(y::Array{T, 1}, w::Array{T2, 1}; halfwin=1, d=2) where {T <: Real, T2 <: Real}
    # constrain the w_min, unless it will lead to matrix division erorr
    w = deepcopy(w) 
    w[w .< 1e-4] .= 1e-4

    n = length(y)
    frame = halfwin*2 + 1;
    if (sum(w) == n); return SG(y; halfwin=halfwin, d=d); end
    
    m = init_mSG(halfwin, d)
    m.S = sgmat_S(halfwin, d); # static

    sgmat_wB!(w[1:frame], m); 
    y_head = @views(m.B[1:halfwin+1, :] * y[1:frame])[:, 1];
    
    y_mid = zeros(T, n - frame - 1) # vector
    @inbounds for i = 1:n-frame-1
        sgmat_wB!(w[i+1:i+frame], m);
        y_mid[i] = @views(dot( m.B[halfwin+1, :], y[i+1:i+frame] ));
    end
    
    sgmat_wB!(w[n-frame+1:n], m);
    y_tail = @views( m.B[halfwin+1:frame, :] * y[n-frame+1:n])[:, 1];
    [y_head; y_mid; y_tail]
end

# weighted Savitzky Golay filter
function wSG_low(y::Array{T, 1}, w::Array{T2, 1}; halfwin=1, d=2) where {T <: Real, T2 <: Real}
    # constrain the w_min, unless it will lead to matrix division erorr
    w = deepcopy(w) 
    w[w .< 1e-4] .= 1e-4

    n = length(y)
    frame = halfwin*2 + 1;
    if (sum(w) == n); return SG(y; halfwin=halfwin, d=d); end

    S = sgmat_S(halfwin, d);
    smat = ones(T2, size(S));
    
    B = sgmat_wB(S, w[1:frame], smat);
    y_head = @views(B[1:halfwin+1, :] * y[1:frame])[:, 1];

    y_mid = zeros(T, n - frame - 1)
    @inbounds for i = 1:n-frame-1
        B = sgmat_wB(S, w[i+1:i+frame], smat);
        y_mid[i] = dot( @view(B[halfwin+1, :]), @view y[i+1:i+frame] );
    end

    B = sgmat_wB(S, w[n-frame+1:n], smat);
    y_tail = @views( B[halfwin+1:frame, :] * y[n-frame+1:n])[:, 1];
    [y_head; y_mid; y_tail]
end

export wSG, wSG_low
