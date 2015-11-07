# Node
abstract AbstractNode{T,N}
abstract AbstractFourierNode{T,N} <: AbstractNode{T,N}

immutable RealFourierNode{T<:FFTW.fftwReal,K,N} <: AbstractFourierNode{T,N}
    data::Array{T,N}
    data_ft::Array{Complex{T},N}
    forward_plan::FFTW.rFFTWPlan{T,K,false,N}
    ranges::NTuple{N,PathRange}
end

immutable ComplexFourierNode{T<:FFTW.fftwComplex,K,N} <:
        AbstractFourierNode{T,N}
    data::Array{T,N}
    data_ft::Array{T,N}
    forward_plan::FFTW.cFFTWPlan{T,K,false,N}
    ranges::NTuple{N,PathRange}
end

immutable InverseFourierNode{T<:FFTW.fftwComplex,R<:FFTW.fftwReal,K,N} <:
        AbstractFourierNode{T,N}
    data::Array{T,N}
    inverse_plan::Base.DFT.ScaledPlan{T,FFTW.cFFTWPlan{T,K,true,N},R}
    ranges::NTuple{N,PathRange}
end

immutable ComplexInverseFourierNode{T<:FFTW.fftwComplex,R<:FFTW.fftwReal,K,N} <:
        AbstractFourierNode{T,N}
    data::Array{T,N}
    data_ft::Array{T,N}
    forward_plan::FFTW.cFFTWPlan{T,K,false,N}
    inverse_plan::Base.DFT.ScaledPlan{T,FFTW.cFFTWPlan{T,K,true,N},R}
    ranges::NTuple{N,PathRange}
end

function AbstractFourierNode{T<:FFTW.fftwReal,N}(
        data::Array{T,N}, fourierdims::Vector{Int}, ranges::NTuple{N,PathRange};
        flags = FFTW.ESTIMATE, timelimit = Inf)
    plan = plan_rfft(data, fourierdims ; flags = flags, timelimit = timelimit)
    data_ft = plan * data
    RealFourierNode(data, data_ft, plan, ranges)
end
function AbstractFourierNode{T<:FFTW.fftwComplex,N}(data::Array{T,N},
        fourierdims::Vector{Int}, ranges::NTuple{N,PathRange};
        flags = FFTW.ESTIMATE, timelimit = Inf)
    plan = plan_fft(data, fourierdims ; flags = flags, timelimit = timelimit)
    data_ft = plan * data
    ComplexFourierNode(data, data_ft, plan, ranges)
end
function AbstractFourierNode{T<:FFTW.fftwNumber,N}(
        data::Array{T,N}, fourierdims::Vector{Int},
        subscripts::NTuple{N, PathKey}; args...)
    ranges =
        ntuple(k -> PathRange(subscripts[k] => (1:1:size(data,k))), ndims(data))
    AbstractFourierNode(data, fourierdims, ranges; args...)
end
AbstractFourierNode(data, fourierdims::Int, subscripts ; args...) =
    AbstractFourierNode(data, collect(fourierdims), subscripts ; args...)

Base.fft!(node::AbstractFourierNode) =
    A_mul_B!(node.data_ft, node.forward_plan, node.data)

Base.ifft!{T<:FFTW.fftwComplex}(node::AbstractFourierNode{T}) =
    A_mul_B!(node.data, node.inverse_plan, node.data)

function pathdepth(node::AbstractFourierNode, refkey::PathKey)
    mapreduce(p -> pathdepth(p.first, refkey), max, 1, node.ranges)
end

function transform!(node_in::AbstractFourierNode, node_out::AbstractFourierNode,
                    ψ::Analytic1DFilter)
    inds = fill!(Array(Union{Colon,Int}, ndims(node_in.data)), Colon())
    N = length(node_in.forward_plan.region[1])
    ψlast = ψ.posfirst + length(ψ.pos) - 1
    # Positive frequencies, excluding midpoint
    @inbounds for ω in ψ.posfirst:max(N>>1 - 1, ψlast)
        inds[node_in.forward_plan.region[1]] = 1 + ω
        view_in = ArrayViews.view(node_in.data_ft, inds...)
        view_out = ArrayViews.view(node_out.data, inds...)
        @inbounds for id in eachindex(view_in)
            view_out[id] = view_in[id] * ψ[1 + ω]
        end
    end
end

function transform!(node_in::RealFourierNode, node_out::AbstractFourierNode,
                    ψ::FullResolution1DFilter)
    inds = fill!(Array(Union{Colon,Int}, ndims(node_in.data)), Colon())
    N = length(node_in.forward_plan.region[1])
    # Positive frequencies, including midpoint
    @inbounds for ω in 1:(N>>1)
        inds[node_in.forward_plan.region[1]] = 1 + ω
        view_in = ArrayViews.view(node_in.data_ft, inds...)
        view_out = ArrayViews.view(node_out.data, inds...)
        @inbounds for id in eachindex(view_in)
            view_out[id] = view_in[id] * ψ[1 + ω]
        end
    end
    # Negative frequencies
    @inbounds for ω in (N>>1+1):(N-1)
        inds[node_in.forward_plan.region[1]] = 1 + N - ω
        view_in = ArrayViews.view(node_in.data_ft, inds...)
        inds[node_in.forward_plan.region[1]] = 1 + ω
        view_out = ArrayViews.view(node_out.data, inds...)
        @inbounds for id in eachindex(view_in)
            view_out[id] = view_in[id] * ψ[1 + ω]
        end
    end
end

function transform!(node_in::ComplexFourierNode, node_out::AbstractFourierNode,
                    ψ::FullResolution1DFilter)
    inds = fill!(Array(Union{Colon,Int}, ndims(node_in.data)), Colon())
    N = length(node_in.forward_plan.region[1])
    # All nonzero frequencies in a single loop
    @inbounds for ω in 1:(N-1)
        inds[node_in.forward_plan.region[1]] = 1 + ω
        view_in = ArrayViews.view(node_in.data_ft, inds...)
        view_out = ArrayViews.view(node_out.data, inds...)
        @inbounds for id in eachindex(view_in)
            view_out[id] = view_in[id] * ψ[1 + ω]
        end
    end
end

function transform!(node_in::RealFourierNode, node_out::AbstractFourierNode,
                    ψ::Vanishing1DFilter)
    inds = fill!(Array(Union{Colon,Int}, ndims(node_in.data)), Colon())
    N = length(node_in.forward_plan.region[1])
    # Positive frequencies, excluding midpoint
    @inbounds for ω in 1:(N>>1)
        inds[node_in.forward_plan.region[1]] = 1 + ω
        view_in = ArrayViews.view(node_in.data_ft, inds...)
        view_out = ArrayViews.view(node_out.data, inds...)
        @inbounds for id in eachindex(view_in)
            view_out[id] = view_in[id] * ψ[1 + ω]
        end
    end
end
