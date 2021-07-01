module phenofit

using StaticArrays: zeros
using Dates: floor
using Plots
# using JLD2

export  whittaker2,
        whittaker2!, 
        smooth2!, 
        smooth2_c!, 
        wTSM, 
        wBisquare
        # whittaker1,
        # whittaker1!,

include("base/main_Ipaper.jl")

# include("raster/raster.jl")
# include("raster/ncread2.jl")
include("QC/qc_FparLai.jl")

# include("smooth_whittaker/whit2_cpp.jl")
# include("smooth_whittaker/whittaker2.jl")
include("smooth_whittaker/lambda_init.jl")
# include("smooth_whittaker/lambda_cv.jl")
include("smooth_whittaker/lambda_vcurve.jl")
include("smooth_whittaker/whit2.jl")
include("smooth_whittaker/smooth_whit.jl")
include("smooth_whittaker/smooth_whit_GEE.jl")
# include("smooth_whittaker/smooth_SG.jl")
include("smooth_SG/smooth_SG.jl")

include("season/season.jl")

include("weights/wBisquare.jl")
include("weights/wTSM.jl")

include("curvefit/plot_input.jl")
include("curvefit/spike.jl")
include("curvefit/movmean.jl")

end
