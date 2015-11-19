# Node
abstract AbstractNode{T<:Number,N}

immutable RealFourierNode{T<:FFTW.fftwReal,N} <: AbstractNode{T,N}
    data::Array{Complex{T},N}
    forwardplan::FFTW.rFFTWPlan{T,-1,false,N}
    ranges::NTuple{N,PathRange}
end

immutable ComplexFourierNode{T<:FFTW.fftwComplex,N} <: AbstractNode{T,N}
    data::Array{T,N}
    forwardplan::FFTW.cFFTWPlan{T,-1,false,N}
    ranges::NTuple{N,PathRange}
end

immutable SpatialNode{T<:Number,N} <: AbstractNode{T,N}
    data::Array{T,N}
    ranges::NTuple{N,PathRange}
end

immutable InverseFourierNode{T<:FFTW.fftwComplex,N,R<:FFTW.fftwReal} <:
        AbstractNode{T,N}
    data::Array{T,N}
    inverseplan::Base.DFT.ScaledPlan{R,FFTW.cFFTWPlan{T,1,false,N},T}
    ranges::NTuple{N,PathRange}
end

function setup{T<:FFTW.fftwReal,N}(
        data::Array{T,N},
        fourierdims::Vector{Int},
        ranges::NTuple{N,PathRange};
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    forwardplan =
        plan_rfft(data, fourierdims ; flags = flags, timelimit = timelimit)
    data_ft = plan * data
    RealFourierNode(data, data_ft, plan, inverseplans, ranges)
end
function setup{T<:FFTW.fftwComplex,N}(
        data::Array{T,N},
        fourierdims::Vector{Int},
        ranges::NTuple{N,PathRange};
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    plan = plan_fft(data, fourierdims ; flags = flags, timelimit = timelimit)
    data_ft = plan * data
    ComplexFourierNode(data, data_ft, plan, ranges)
end
function setup{T<:FFTW.fftwComplex,N}(
        data::Array{T,N},
        fourierdims::Vector{Int},
        subscripts::NTuple{N,PathKey};
        args...)
    ranges =
        ntuple(k -> PathRange(subscripts[k] => (1:1:size(data,k))), ndims(data))
    setup(data, fourierdims, ranges; args...)
end
setup(data, fourierdims::Int, subscripts ; args...) =
    setup(data, collect(fourierdims), subscripts ; args...)

function pathdepth(node::AbstractNode, refkey::PathKey)
    mapreduce(p -> pathdepth(p.first, refkey), max, 1, node.ranges)
end

function transform!(
        node_in::AbstractNode{FourierDomain{1}},
        node_out::AbstractNode{FourierDomain{1}},
        ψ::Analytic1DFilter)
    inds = fill!(Array(Union{Colon,Int}, ndims(node_in.data)), Colon())
    N = length(node_in.forward_plan.region[1])
    ψlast = ψ.posfirst + length(ψ.pos) - 1
    # Positive frequencies, excluding midpoint
    @inbounds for ω in ψ.posfirst:max(N>>1 - 1, ψlast)
        inds[node_in.forwardplan.region[1]] = 1 + ω
        view_in = ArrayViews.view(node_in.data_ft, inds...)
        view_out = ArrayViews.view(node_out.data, inds...)
        @inbounds for id in eachindex(view_in)
            view_out[id] = view_in[id] * ψ[1 + ω]
        end
    end
end

function transform!(
        node_in::RealFourierNode{FourierDomain{1}},
        node_out::AbstractNode{FourierDomain{1}},
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

function transform!(
        node_in::ComplexFourierNode,
        node_out::AbstractNode{FourierDomain{1}},
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

function transform!(
        node_in::RealFourierNode,
        node_out::AbstractFourierNode{Fourier},
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
