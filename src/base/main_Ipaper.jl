using Glob
using Dates

# x and y is equal length at here
# - `I_y`: which y is corresponding to current x[i]
function match2(x, y)
    I_y = ones(Int, length(x))
    for i in 1:length(x)
        I_y[i] = findfirst(x[i] .== y)
        # I_x[i] = findfirst(y[i] .== x)
    end
    I_y # 
end

# function match2(x, y)
#     I_x = ones(Int, length(x))
#     I_y = ones(Int, length(y))
#     for i in 1:length(x)
#         I_x[i] = findfirst(x[i] .== y)
#     end
#     for i in 1:length(y)
#         I_y[i] = findfirst(y[i] .== x)
#     end
#     I_x, I_y
# end
function list_dir(indir)
    filter(x -> isdir(joinpath(indir, x)), readdir(indir, join = true))
end


function str_extract(x, pattern = r"\d{4}_\d{1}_\d") 
    if typeof(x) == String; x = [x]; end
    str = match.(pattern, x)
    str = map(x -> x.match, str)
end


function split(list, names) 
    grps = unique(names)
    map(grp -> Pair(grp, list[names .== grp]), grps)
end


function CartesianIndex2Int(x, ind)
    I = LinearIndices(x)
    # I = 1:prod(size(x))
    I[ind]
end


"""
    merge_pdf("*.pdf", output="Plot.pdf")

merge multiple pdf files by `pdftk`
"""
function merge_pdf(input, output="Plot.pdf"; is_del = false)
    # input = abspath(input)
    files = glob(input)
    id = str_extract(basename.(files), r"\d{1,}") 
    id = parse.(Int32, id) |> sortperm
    files = files[id]

    run(`pdftk $files cat output $output`)
    if is_del; run(`rm $files`); end
    nothing
end


"""
    open pdf file in SumatraPDF
"""
macro show_pdf(file)
    run(`/mnt/c/WINDOWS/SumatraPDF.exe $file`; wait = false)
    nothing
end

macro show_file(file)
    run(`cmd /c $file`; wait = false)
    nothing
end

macro methods(func)
    :(methods($func))
end

macro savefig(plt::Plots.Plot, fn::AbstractString, show = true)
    :(savefig($plt, $fn))
    if (show); :(@show_file($fn)); end
end

macro savefig(fn::AbstractString, show = true)
    :(savefig($fn))
    println("jldsave")
    if (show); :(@show_file($fn)); end
end

# # import PyCall
# # PyPlot
# # pdf = PyCall.pyimport("matplotlib.backends.backend_pdf")
# function write_pdf(figures, file = "Plot.pdf")
#     if isfile(file); rm(file); end
#     # pdf()
#     pdffile = pdf.PdfPages(file) # create pdf file
#     [pdffile.savefig(f) for f in figures] # add figures to file
#     pdffile.close() # close pdf file    
# end
get_dn(date, days = 8) = fld.(Dates.dayofyear.(date) .- 1, days) .+ 1

set_value!(x, con, value) = begin
    x[con] .= value
    Nothing
end

export match2, split, str_extract, list_dir, CartesianIndex2Int, merge_pdf, 
    @show_file, @show_pdf, 
    @methods,
    @savefig,
    get_dn, 
    set_value!
