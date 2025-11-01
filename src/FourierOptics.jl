module FourierOptics

import FromFile: @from
@from "Modules/Bare.jl" import Bare
@from "Modules/Propagation.jl" import Propagation

using Reexport
@reexport using .Bare
@reexport using .Propagation
@reexport using .Fourier

export Bare, Propagation, Fourier

end
