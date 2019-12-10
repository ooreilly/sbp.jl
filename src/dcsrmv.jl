module DCSRMV
"""
Export matrices and vectors to disk so that they can be read by the code dcsrmv.
This code performs sparse CSR matrix- dense vector multiplication on the GPU.

"""

using SparseArrays

vector_type = [Int32, Float32, Float64]
enum_Int32 = Int32(1)
enum_Float32 = Int32(2)
enum_Float64 = Int32(3)


function write_vector(filename::String, u::Array{Int32, 1}; verbose::Bool=false)
        write_vector(filename, u, enum_Int32, verbose=verbose)
end

function write_vector(filename::String, u::Array{Float32, 1};
                                        verbose::Bool=false)
       write_vector(filename, u, enum_Float32, verbose=verbose)
end

function write_vector(filename::String, u::Array{Float64, 1};
                                        verbose::Bool=false)
       write_vector(filename, u, enum_Float64, verbose=verbose)
end

function write_vector(filename::String, u::Array{T, 1}, enum::Int32;
                                        verbose::Bool=false) where T
        fh = open(filename, "w")
        write(fh, enum)
        n = length(u)
        write(fh, n)
        write(fh, u)
        close(fh)
        if verbose
                println("Wrote vector of length ", n, " to ", filename)
        end
end

function read_vector(filename::String; verbose::Bool=false)
        fh = open(filename, "r")
        enum = read(fh, Int32)
        n = read(fh, Int64)
        type = vector_type[enum]
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

function write_csr_matrix(filename::String, A::SparseMatrixCSC{Float64, Int64}; 
                          verbose::Bool=false)
        fh = open(filename, "w")
        B = copy(A')
        write(fh, B.m)
        write(fh, B.n)
        nnz = length(B.nzval)
        write(fh, nnz)

        # Convert to zero-indexing
        rows =  B.colptr .- 1
        cols =  B.rowval .- 1
        write(fh, rows)
        write(fh, cols)
        write(fh, B.nzval)
        close(fh)

        if verbose
                println("Wrote CSR matrix <", B.m, " x ", B.n, "> containing ",
                        nnz, " non-zeros to ", filename) 
        end
end

function read_csr_matrix(filename::String; verbose=false::Bool)
        fh = open(filename, "r")
        m = read(fh, Int64)
        n = read(fh, Int64)
        nnz = read(fh, Int64)

        colptr = zeros(Int64, m + 1)
        for i=1:(m + 1)
                colptr[i] = read(fh, Int64) + 1
        end

        rowval = zeros(Int64, nnz)
        for i=1:nnz
                rowval[i] = read(fh, Int64) + 1
        end

        nzval = zeros(Float64, nnz)
        for i=1:nnz
                nzval[i] = read(fh, Float64)
        end
        close(fh)

        B = SparseMatrixCSC(m, n, colptr, rowval, nzval)
        A = copy(B')

        if verbose
                println("Read CSR matrix <", B.m, " x ", B.n, "> containing ",
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
                        A[i,j] = read(fh, Float64)
                end
        end

        if verbose
                println("Read dense matrix <", m, " x ", n, "> to ", filename) 
        end

        return A

end

end
