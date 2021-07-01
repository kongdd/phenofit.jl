using DataFrames

is_empty(x) = length(x) == 0

# throw(ArgumentError("syntax df[column] is not supported use df[!, column] instead"))
Base.getindex(d::AbstractDataFrame, col::Union{Symbol, Integer, AbstractString}) = begin
    # println("hello")
    d[:, col]
end

Range(x::AbstractArray) = maximum(x) - minimum(x)


include("findpeaks.jl")
include("findpeaks_season.jl")
include("check_season.jl")

export findpeaks, findpeaks_season, meltPeakTrough, 
    is_empty, Range
