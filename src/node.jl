# Node
abstract AbstractNode{T,N}
abstract AbstractFourierNode{T,N} <: AbstractNode{T,N}

Base.complex{T<:Real}(::Type{T}) = Complex{T}
Base.complex{T<:Complex}(::Type{T}) = T

immutable RealFourierNode{T<:Base.FFTW.fftwReal,K,N} <: AbstractFourierNode{T,N}
    data::Array{T,N}
    data_ft::Array{Complex{T},N}
    forward_plan::Base.DFT.FFTW.rFFTWPlan{T,K,false,N}
    ranges::NTuple{N,PathRange}
end

immutable ComplexFourierNode{T<:Base.FFTW.fftwReal,N} <:
        AbstractFourierNode{Complex{T},N}
    data::Array{Complex{T},N}
    data_ft::Array{Complex{T},N}
    forward_plan::Base.DFT.FFTW.cFFTWPlan{Complex{T},K,false,N}
    ranges::NTuple{N,PathRange}
end

immutable InverseFourierNode{T<:Base.FFTW.fftwReal,K,N} <:
        AbstractFourierNode{Complex{T},N}
    data::Array{Complex{T},N}
    inverse_plan::Base.DFT.ScaledPlan{Complex{T},
        Base.DFT.FFTW.cFFTWPlan{Complex{T},K,true,N},T}
    ranges::NTuple{N,PathRange}
end

immutable ComplexInverseFourierNode{T<:Base.FFTW.fftwReal,K,N} <:
        AbstractFourierNode{Complex{T},N}
    data::Array{Complex{T},N}
    data_ft::Array{Complex{T},N}
    forward_plan::Base.DFT.FFTW.cFFTWPlan{Complex{T},K,false,N}
    inverse_plan::Base.DFT.ScaledPlan{Complex{T},
        Base.DFT.FFTW.cFFTWPlan{Complex{T},K,true,N},T}
    ranges::NTuple{N,PathRange}
end

function AbstractFourierNode{T<:Real,N}(data::Array{T,N},
                                        fourierdims::Vector{Int},
                                        ranges::NTuple{N,PathRange};
                                        flags = FFTW.ESTIMATE,
                                        timelimit = Inf)
    plan = plan_rfft(data, fourierdims ; flags = flags, timelimit = timelimit)
    data_ft = plan * data
    RealFourierNode(data, data_ft, plan, ranges)
end
function AbstractFourierNode{T<:Complex,N}(data::Array{T,N},
                                           fourierdims::Vector{Int},
                                           ranges::NTuple{N,PathRange};
                                           flags = FFTW.ESTIMATE,
                                           timelimit = Inf)
    plan = plan_fft(data, fourierdims ; flags = flags, timelimit = timelimit)
    data_ft = plan * data
    ComplexFourierNode{T,N}(data, data_ft, plan, ranges)
end
function AbstractFourierNode{T<:Number,N}(data::Array{T,N},
                                          fourierdims::Vector{Int},
                                          subscripts::NTuple{N, PathKey};
                                          args...)
    ranges =
        ntuple(k -> PathRange(subscripts[k] => (1:1:size(data,k))), ndims(data))
    AbstractFourierNode(data, fourierdims, ranges; args...)
end
AbstractFourierNode(data, fourierdims::Int, subscripts ; args...) =
    AbstractFourierNode(data, collect(fourierdims), subscripts ; args...)

Base.fft!(node::AbstractFourierNode) =
    A_mul_B!(node.data_ft, node.forward_plan, node.data)

Base.ifft!{T<:Complex}(node::AbstractFourierNode{T}) =
    A_mul_B!(node.data, node.inverse_plan, node.data)
