
diffs(x::Array{T,1}, d = 2) where T <: Real = d == 0 ? x : diffs(diff(x), d-1) 

fidelity(y, z, w) = log(sum(@.(w*(y-z)^2)))
roughness(z, d = 2) = log(sum( diffs(z, d) .^ 2 ))

"""
    lambda_vcurve(y::Array{T,1}, w::Array{T2,1}; 
    is_plot = false, lg_lambda_min = 0.1, lg_lambda_max = 3)

# Return
- `lambda`: optimal lambda
"""
function lambda_vcurve(y::Array{T,1}, w::Array{T2,1}; 
    is_plot = false, lg_lambda_min = 0.1, lg_lambda_max = 3) where { 
        T <: Real, T2 <: Real }

    lg_lambdas = lg_lambda_min:0.1:lg_lambda_max
    n = length(lg_lambdas)
    fits = zeros(n)
    pens = zeros(n)
    
    # least of memory used
    m = length(y)
    z = zeros(T, m)
    c = zeros(Float32, m)
    d = zeros(Float32, m)
    e = zeros(Float32, m)

    for i in 1:n
        lambda = 10^lg_lambdas[i]
        # z, cve = whit2(y, w, lambda, include_cve = false)
        cve = whit2!(y, w, lambda, z, c, d, e, include_cve = false)
        # println(length(y), length(z), length(w))
        fits[i] = fidelity(y, z, w)
        pens[i] = roughness(z)
    end
    
    dfits = diff(fits)
    dpens = diff(pens)
    
    llastep = lg_lambdas[2] - lg_lambdas[1] 
    v = @. sqrt(dfits^2 + dpens^2)/(log(10) * llastep)
    lamids = (lg_lambdas[2:end] + lg_lambdas[1:end-1])/2
    k = argmin(v)
    opt_lambda = 10^lamids[k]
    # z = whit2(y, lambda, w)
    if is_plot; plot_lambda(y, w, lamids, v) |> display; end
    opt_lambda
end

"""
    lambda_cv(y::Array{T,1}, w::Array{T2,1}; 
    is_plot = false, lg_lambda_min = 0.1, lg_lambda_max = 3)

# Return
- `lambda`: optimal lambda
"""
function lambda_cv(y::Array{T,1}, w::Array{T2,1}; 
    is_plot = false, lg_lambda_min = 0.1, lg_lambda_max = 3) where { 
        T <: Real, T2 <: Real }

    lg_lambdas = lg_lambda_min:0.1:lg_lambda_max
    n = length(lg_lambdas)

    # least of memory used
    m = length(y)
    z = zeros(T, m)
    c = zeros(Float32, m)
    d = zeros(Float32, m)
    e = zeros(Float32, m)
    cvs = zeros(n)
    # pens = zeros(n)
    for i in 1:n
        lambda = 10^lg_lambdas[i]
        # z, cvs[i] = whit2(y, w, lambda)
        cvs[i] = whit2!(y, w, lambda, z, c, d, e, include_cve = true)
        # fits[i] = fidelity(y, z, w)
        # pens[i] = roughness(z, d)
    end
    # dfits = diff(fits)
    # dpens = diff(pens)
    # llastep = lg_lambdas[2] - lg_lambdas[1] 
    # v = @. sqrt(dfits^2 + dpens^2)/(log(10) * llastep)
    # lamids = (lg_lambdas[2:end] + lg_lambdas[1:end-1])/2
    k = argmin(cvs)
    opt_lambda = 10^lg_lambdas[k]
    if is_plot; plot_lambda(y, w, lg_lambdas, cvs) |> display; end
    opt_lambda
end

# x: lambdas candidates
function plot_lambda(y, w, lg_lambdas, cvs)
    k = argmin(cvs)
    opt_lambda = 10^lg_lambdas[k]

    p_v = plot(lg_lambdas, cvs, label = "Generalized CV", frame = :box)
    scatter!(p_v, lg_lambdas, cvs, legend = false)
    scatter!(p_v, [lg_lambdas[k]], [cvs[k]], 
        m = (10, :transparent, stroke(1, "red")),
        legend = false)
    vline!(p_v, [lg_lambdas[k]], color = "red", linestyle = :dash)
    
    xlim = (0, length(y))
    p1   = plot(y, xlim = xlim, frame = :box)
    z_2,  = whit2(y, w, 2.0)
    z_15, = whit2(y, w, 15.0)
    z_opt, = whit2(y, w, opt_lambda)
    plot!(p1, z_2, label = "lambda = 2")
    plot!(p1, z_15, label = "lambda = 15")
    plot!(p1, z_opt, label = "lambda = $(round(opt_lambda, digits = 3))")
    plot(p_v, p1, layout = (1, 2), size = (700, 480))
end

export lambda_vcurve, lambda_cv, plot_lambda
