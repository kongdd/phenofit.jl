using Pkg
pkg"activate .."
using Documenter
using phenofit

makedocs(sitename="My Documentation", 
    modules = [phenofit])
