# Node
abstract AbstractNode{T<:Number,N}
abstract AbstractFourierNode{T<:Number,N} <: AbstractNode{T,N}

AbstractFourierNode{T<:FFTW.fftwReal}(node::AbstractNode{T}, args...) =
    RealFourierNode(node, args...)
AbstractFourierNode{T<:FFTW.fftwComplex}(node::AbstractNode{T}, args...) =
    ComplexFourierNode(node, args...)

immutable RealFourierNode{T<:FFTW.fftwComplex,R<:FFTW.fftwReal,N} <:
        AbstractFourierNode{T,N}
    data::Array{T,N}
    forwardplan::FFTW.rFFTWPlan{R,-1,false,N}
    ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}}
end

function RealFourierNode{T<:FFTW.fftwReal}(
        node::AbstractNode{T},
        region::AbstractArray{Int,1},
        flags::UInt32,
        timelimit = Inf)
    forwardplan = plan_rfft(
            node.data, region, flags = flags, timelimit = timelimit)
    RealFourierNode(forwardplan * node.data, forwardplan, node.ranges)
end

immutable ComplexFourierNode{T<:FFTW.fftwComplex,N} <:
        AbstractFourierNode{T,N}
    data::Array{T,N}
    forwardplan::FFTW.cFFTWPlan{T,-1,false,N}
    ranges::NTuple{N,Pair{PathKey,StepRange{Int,Int}}}
end

function ComplexFourierNode{T<:FFTW.fftwComplex}(
        node::AbstractNode{T},
        region::AbstractArray{Int,1};
        flags = FFTW.ESTIMATE,
        timelimit = Inf)
    forwardplan =
        plan_fft(node.data, region ; flags = flags, timelimit = timelimit)
    ComplexFourierNode(forwardplan * node.data, forwardplan, node.ranges)
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
    InvComplexFourierNode(inverseplan * node.data, inverseplan, node.ranges)
end

function pathdepth(node::AbstractNode, refkey::PathKey)
    mapreduce(p -> pathdepth(p.first, refkey), max, 1, node.ranges)
end
