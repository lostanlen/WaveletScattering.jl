# Node
abstract AbstractNode{T<:Number,N}
abstract AbstractFourierNode{T<:Number,N} <: AbstractNode{T,N}

immutable RealFourierNode{
        T<:FFTW.fftwComplex,R<:FFTW.fftwReal,N} <: AbstractFourierNode{T,N}
    data::Array{T,N}
    forwardplan::FFTW.rFFTWPlan{R,-1,false,N}
    ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}}
end

function RealFourierNode{T<:FFTW.fftwReal,N}(
        data::Array{T,N},
        region::Vector{Int},
        pathkeys::NTuple{N,PathKey};
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    ranges = ntuple(k -> (pathkeys[k] => (1:1:size(data,k))), ndims(data))
    forwardplan =
        plan_rfft(data, region ; flags = flags, timelimit = timelimit)
    RealFourierNode(forwardplan * data, forwardplan, ranges)
end

function RealFourierNode{T<:FFTW.fftwComplex}(
        node::AbstractNode{T},
        region::Vector{Int};
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    forwardplan =
        plan_rfft(data, region ; flags = flags, timelimit = timelimit)
    RealFourierNode(forwardplan * data, forwardplan, node.ranges)
end

immutable ComplexFourierNode{T<:FFTW.fftwComplex,N} <: AbstractFourierNode{T,N}
    data::Array{T,N}
    forwardplan::FFTW.cFFTWPlan{T,-1,false,N}
    ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}}
end

function ComplexFourierNode{T<:FFTW.fftwComplex,N}(
        data::Array{T,N},
        region::Vector{Int},
        subscripts::NTuple{N,PathKey};
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    ranges = ntuple(k -> (subscripts[k] => (1:1:size(data,k))), ndims(data))
    forwardplan =
        plan_fft(data, region ; flags = flags, timelimit = timelimit)
    ComplexFourierNode(forwardplan * data, forwardplan, ranges)
end

function ComplexFourierNode{T<:FFTW.fftwComplex}(
        node::AbstractNode{T},
        region::Vector{Int};
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    forwardplan =
        plan_fft(data, region ; flags = flags, timelimit = timelimit)
    ComplexFourierNode(forwardplan * data, forwardplan, node.ranges)
end

immutable Node{T<:Number,N} <: AbstractNode{T,N}
    data::Array{T,N}
    ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}}
end

immutable InvComplexFourierNode{T<:FFTW.fftwComplex,N,R<:FFTW.fftwReal} <:
        AbstractNode{T,N}
    data::Array{T,N}
    inverseplan::Base.DFT.ScaledPlan{R,FFTW.cFFTWPlan{T,1,false,N},T}
    ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}}
end

function InvComplexFourierNode{T<:FFTW.fftwComplex}(
        node::AbstractNode{T},
        region::Vector{Int} ;
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    inverseplan =
        plan_ifft(node.data, region ; flags = flags, timelimit = timelimit)
    InvComplexFourierNode(inverseplan * data, inverseplan, node.ranges)
end

function pathdepth(node::AbstractNode, refkey::PathKey)
    mapreduce(p -> pathdepth(p.first, refkey), max, 1, node.ranges)
end

function transform!(
        node_out::AbstractNode,
        node_in::ComplexFourierNode,
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
