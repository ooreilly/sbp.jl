module lmv
"""
Export matrices and vectors to disk so that they can be read by the code dcsrmv.
This code performs sparse CSR matrix- dense vector multiplication on the GPU.

"""

using SparseArrays

function write_vector(filename::String, u::Array{T, 1};
                                        verbose::Bool=false) where T
        fh = open(filename, "w")
        write(fh, sizeof(T))
        n = length(u)
        bytes = write(fh, UInt64(n))
        write(fh, u)
        close(fh)
        if verbose
                println("Wrote vector of length ", n, " to ", filename)
        end
end

function read_vector(type, filename::String; verbose::Bool=false)
        fh = open(filename, "r")
        type_size = read(fh, UInt64)
        n = read(fh, UInt64)
        @assert(sizeof(type) == type_size)
        u = zeros(type, n)
        for i=1:n
                @inbounds u[i] = read(fh, type)
        end
        close(fh)
        if verbose
                println("read vector of length ", n, " from ", filename)
        end
        return u
end

function write_csr_matrix(filename::String, A::SparseMatrixCSC{Tv, Ti}; 
                          verbose::Bool=false) where {Tv, Ti}
        fh = open(filename, "w")
        B = copy(A')
        write(fh, sizeof(Tv))
        write(fh, sizeof(Ti))
        write(fh, UInt64(B.n))
        write(fh, UInt64(B.m))
        nnz = length(B.nzval)
        write(fh, UInt64(nnz))

        # Convert to zero-indexing
        rows =  B.colptr .- 1
        cols =  B.rowval .- 1

        write(fh, Ti.(rows))
        write(fh, Ti.(cols))
        write(fh, Tv.(B.nzval))
        close(fh)

        if verbose
                println("Wrote CSR matrix <", A.m, " x ", A.n, "> containing ",
                        nnz, " non-zeros to ", filename) 
        end
end

function read_csr_matrix(Tv, Ti, filename::String; verbose=false::Bool)
        fh = open(filename, "r")
        tv = read(fh, UInt64)
        ti = read(fh, UInt64)
        m = read(fh, UInt64)
        n = read(fh, UInt64)
        nnz = read(fh, UInt64)

        colptr = zeros(Ti, m + 1)
        for i=1:(m + 1)
                colptr[i] = read(fh, Ti) + 1
        end

        rowval = zeros(Ti, nnz)
        for i=1:nnz
                rowval[i] = read(fh, Ti) + 1
        end

        nzval = zeros(Tv, nnz)
        for i=1:nnz
                nzval[i] = read(fh, Tv)
        end
        close(fh)

        B = SparseMatrixCSC(n, m, colptr, rowval, nzval)
        A = copy(B')

        if verbose
                println("Read CSR matrix <", A.m, " x ", A.n, "> containing ",
                        nnz, " non-zeros from ", filename) 
        end

        return A
end

function write_dense_matrix(filename::String, A::Array{Float64, 2};
                                              verbose=false::Bool)

        fh = open(filename, "w")
        m = size(A, 1)
        n = size(A, 2)
        write(fh, m)
        write(fh, n)
        write(fh, A)
        close(fh)

        if verbose
                println("Wrote dense matrix <", m, " x ", n, "> to ", filename) 
        end

end

function read_dense_matrix(filename::String; verbose=false::Bool)

        fh = open(filename, "r")
        m = read(fh, Int64)
        n = read(fh, Int64)
        A = zeros(m, n)
        for j=1:n
                for i=1:m
                        @inbounds A[i,j] = read(fh, Float64)
                end
        end

        if verbose
                println("Read dense matrix <", m, " x ", n, "> to ", filename) 
        end

        return A

end

function write_config(filename::String, args::Dict; verbose=false::Bool,
                                        skip_types=false::Bool)
        fh = open(filename, "w")
        for arg in args
                if skip_types == true
                        println(fh, arg.first, "=", arg.second)
                else
                        println(fh, arg.first, ":", string(typeof(arg.second)),
                                    "=", arg.second)
                end
        end
        close(fh)

        if verbose
                println("Wrote configuration file: ", filename)
        end
end

function read_config(filename::String; verbose=false::Bool)
        fh = open(filename, "r")
        cfg = Dict()
        for ln in eachline(fh)
                name, tmp = split(ln, ":")
                typ, value = split(tmp, "=")
                typ = eval(Meta.parse(typ))
                cfg[name] = parse(typ, value)
        end
        close(fh)

        if verbose
                println("Read configuration file: ", filename)
        end

        return cfg
end


end
