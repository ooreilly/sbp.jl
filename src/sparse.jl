module Sparse
using SparseArrays


"""
 Allocates a block sparse matrix that stores nnz non-zero entries. 
 The matrix contains only zeros after allocation.

 Input:
    rows: A vector describing the number of elements in each block row the
          block matrix has.
 columns: A vector describing the number of elements in each block column the
          block matrix has.

 Output:
       Z: Sparse matrix (all zeros)

 Example:  
 >> Z = block_matrix([2, 3],[2, 3])
 >> size(Z)
 
 ans =
 
      5     5
"""
function block_matrix(rows::AbstractArray, columns::AbstractArray)
        m = sum(rows);
        n = sum(columns);
        
        Z = spzeros(m,n);
        return Z
end


"""
 Z = block_matrix_insert(B,rows,columns,block)
 Inserts a submatrix into a block matrix at block position (i,j)

 Input:
       B: Block matrix (see block_matrix.m for a description).
    rows: Number of elements in each block row of the block matrix.
 columns: Number of elements in each block row of the block matrix.
     i,j: Block indices denoting where the submatrix should be placed.
   block: The submatrix to insert into the block matrix. 
          The size of the submatrix must match the number given in rows(i)
          and columns(j).

 Example:
 >> Z = block_matrix([2, 3],[2, 3])
 >> B = block_matrix_insert(Z,[2, 3],[2, 3],1,2,ones(2,3))
 >> B
 
 B =
 
    (1,3)        1
    (2,3)        1
    (1,4)        1
    (2,4)        1
    (1,5)        1
    (2,5)        1
>> full(B)

ans =

     0     0     1     1     1
     0     0     1     1     1
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
"""
function block_matrix_insert(B::AbstractArray,
                             rows::AbstractArray,
                             columns::AbstractArray,
                             i::Int64,
                             j::Int64,
                             block::AbstractArray)



@assert(i <= length(rows),          "Index out of bounds.")
@assert(j <= length(columns),       "Index out of bounds.")
@assert(size(block,1) == rows[i],   "Dimension mismatch")
@assert(size(block,2) == columns[j],"Dimension mismatch")

        r = vcat([1],rows)
        c = vcat([1],columns)
        row_offset = sum(r[1:i])
        col_offset = sum(c[1:j])
        B[row_offset:row_offset+r[i+1]-1,col_offset:col_offset+c[j+1]-1] = block
        return B
end


function block_matrix_2x2(a11::AbstractArray, 
                          a12::AbstractArray,
                          a21::AbstractArray,
                          a22::AbstractArray)
        m1 = size(a11, 1)
        n1 = size(a11, 2)
        m2 = size(a22, 1)
        n2 = size(a22, 2)
        @assert size(a12, 1) == size(a11, 1)
        @assert size(a12, 2) == size(a22, 2)
        @assert size(a21, 1) == size(a22, 1)
        @assert size(a21, 2) == size(a11, 2)
        rows = [m1, m2]
        cols = [n1, n2]
        A = block_matrix(rows, cols)
        A = block_matrix_insert(A, rows, cols, 1, 1, a11)
        A = block_matrix_insert(A, rows, cols, 1, 2, a12)
        A = block_matrix_insert(A, rows, cols, 2, 1, a21)
        A = block_matrix_insert(A, rows, cols, 2, 2, a22)
        return A
end

function diag_block_matrix_2x2(a11::AbstractArray,
                               a22::AbstractArray)
        m1 = size(a11, 1)
        n1 = size(a11, 2)
        m2 = size(a22, 1)
        n2 = size(a22, 2)
        rows = [m1, m2]
        cols = [n1, n2]
        A = block_matrix(rows, cols)
        A = block_matrix_insert(A, rows, cols, 1, 1, a11)
        A = block_matrix_insert(A, rows, cols, 2, 2, a22)
        return A
end

end
