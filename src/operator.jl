module Operator
include("config.jl")
using .Config: Tv

using Printf: @printf
using SparseArrays
using DelimitedFiles

major=1
minor=0
patch=0

struct Version
major::Int32
minor::Int32
patch::Int32
end

struct OperatorData
 version::Version
 left::AbstractArray
 right::AbstractArray
 interior::AbstractArray
 offset::Int64
end

struct Operators2D
        nx::Int64
        ny::Int64
        Dx::AbstractArray
        Dy::AbstractArray
        Px::AbstractArray
        Py::AbstractArray
        Hx::AbstractArray
        Hy::AbstractArray
        Bx::AbstractArray
        By::AbstractArray
end

version = Version(1, 0, 0)

function to_matrix(op::OperatorData, m::Int64, n::Int64, is_sparse::Bool=false)
        min_size = size(op.left, 1) + size(op.right, 1)
        if m < min_size
                error("Too few rows")
        end

        if is_sparse
                A = spzeros(m, n)
        else
                A = zeros(m, n)
        end

        # Add left boundary block
        ml = size(op.left, 1)
        nl = size(op.left, 2)
        A[1:ml, 1:nl] = op.left

        # Add right boundary block
        mr = size(op.right, 1)
        nr = size(op.right, 2)
        A[end-mr+1:end, end-nr+1:end] = op.right

        # Add interior stencil
        sm = size(op.interior, 1)
        for i in ml + 1:m - mr
                for j in 1:sm
                        if i + op.offset + j > n
                                break
                        end
                        @inbounds A[i,i+op.offset+j] = op.interior[j]
                end
        end



        return A
end

function to_vector(op::OperatorData, n::Int64)
        min_size = size(op.left, 2) + size(op.right, 2)
        if n < min_size
                error("Too few elements")
        end

        u = zeros(n)

        # Left
        nl = size(op.left, 2)
        u[1:nl] = op.left[1,:]

        # Right
        nr = size(op.right, 2)
        u[end-nr+1:end] = op.right[1,:]

        # Interior
        for i in nl + 1:n - nr
                u[i] = op.interior[1,1]
        end

        return u


end

function parse_version(line::String)
        version = split(line, ".")
        v = version[1][1]
        if v != 'v'
                error("Unknown file format")
        end

        fmajor = parse(Int32, version[1][2:end])
        fminor = parse(Int32, version[2])
        fpatch = parse(Int32, version[3])

        version = Version(fmajor, fminor, fpatch)

        return version
end

function check_version(v::Version)
        if v.major != major || v.minor > minor || v.patch > patch
                error("Incompatible version")
        end
end

function write_operator(filename::String, op::OperatorData)
        f = open(filename, "w")
        @printf(f, "v%d.%d.%d\n\n", 
                op.version.major, op.version.minor, op.version.patch)
        write_array(f, "left", op.left)
        write_array(f, "right", op.right)
        write_stencil(f, "interior", op.interior, op.offset)

        close(f)

end

function read_operator(filename::String)
        f = open(filename)
        line = readline(f)
        version = parse_version(line)
        line = readline(f)
        label, left = read_array(f)
        @assert label == "left"
        label, right = read_array(f)
        @assert label == "right"
        label, interior, offset = read_stencil(f)
        @assert label == "interior"
        close(f)

        return OperatorData(version, left, right, interior, offset)

end

function equals(a::OperatorData, b::OperatorData)
        same = true
        same |= a.version == b.version
        same |= a.offset == b.offset
        same |= isapprox(a.left, b.left)
        same |= isapprox(a.right, b.right)
        same |= isapprox(a.interior, b.interior)
        return same
end

function read_array(f::IOStream)

        label = readline(f)
        line = readline(f)
        sizes = split(line)
        if length(sizes) != 2
                error("Unable to parse array")
        end
        m = parse(Int64, sizes[1])
        n = parse(Int64, sizes[2])
        
        v = zeros(m, n)
        for i in 1:m
                row = readline(f)
                val = split(row)
                for j in 1:n
                        v[i,j] = parse(Tv, val[j])
                end
        end
        line = readline(f)

        return label, v

end

function read_stencil(f::IOStream)

        label = readline(f)
        line = readline(f)
        sizes = split(line)
        if length(sizes) != 2
                error("Unable to parse stencil")
        end
        m = parse(Int64, sizes[1])
        offset = parse(Int64, sizes[2])
        
        v = zeros(m)
        row = readline(f)
        val = split(row)
        for i in 1:m
                v[i] = parse(Tv, val[i])
        end

        line = readline(f)

        return label, v, offset
end

function write_array(f, label, a)
        m = size(a,1)
        n = size(a,2)
        println(f, label)
        println(f, m, " ", n)
        for i in 1:m
                for j in 1:n
                        print(f, a[i,j], " ")
                end
                println(f, "")
        end
        println(f, "")
end

function write_stencil(f, label, a, offset)
        m = size(a,1)
        println(f, label)
        println(f, m, " ", offset)
        for i in 1:m
                print(f, a[i], " ")
        end
        println(f, "")
        println(f, "")
end

end
