module VTK

using Printf: @printf


function vtk_write(filename::String, x::AbstractArray, y::AbstractArray,
                   field::AbstractArray, nx::Int64, ny::Int64)
        numpts = nx * ny
        @assert numpts == size(x, 1)
        @assert numpts == size(y, 1)
        @assert numpts == size(field, 1)
        
        fh = open(filename, "w")

        # Header
        @printf(fh, "# vtk DataFile Version 4.2\n")
        @printf(fh, "vtk output\n")
        @printf(fh, "ASCII\n")
        @printf(fh, "DATASET STRUCTURED_GRID\n")
        @printf(fh, "DIMENSIONS %d %d %d\n", nx, ny, 1)
        @printf(fh, "POINTS %d float\n", numpts)
        
        # Coordinates
        for j in 1:ny
                for i in 1:nx
                        @printf(fh, "%f %f %f\n", 
                                x[j + (i - 1) * ny], 
                                y[j + (i - 1) * ny], 
                                0.0)
                end
        end

        # Field
        @printf(fh, "POINT_DATA %d \n", numpts)
        @printf(fh, "FIELD scalar 1\n")
        @printf(fh, "data 1 %d float\n", numpts)

        for j in 1:ny
                for i in 1:nx
                        @printf(fh, "%f \n", field[j + (i - 1) * ny])
                end
        end

        close(fh)
end

function vtk_read(filename)

        fh = open(filename, "r")

        # Header
        line = readline(fh)
        @assert line == "# vtk DataFile Version 4.2"
        line = readline(fh)
        @assert line == "vtk output"
        line = readline(fh)
        @assert line == "ASCII"
        line = readline(fh)
        @assert line == "DATASET STRUCTURED_GRID"
        line = readline(fh)
        fields = split(line)
        
        # Grid dimensions
        @assert fields[1] == "DIMENSIONS"
        nx = parse(Int64, fields[2])
        ny = parse(Int64, fields[3])

        # Number of points
        line = readline(fh)
        fields = split(line)
        @assert fields[1] == "POINTS"
        numpts = parse(Int64, fields[2])
        @assert fields[3] == "float"

        # Data section
        x = zeros(numpts)
        y = zeros(numpts)
        field = zeros(numpts)

        for j in 1:ny
                for i in 1:nx
                        line = readline(fh)
                        fields = split(line)
                        idx = j + (i - 1) * ny
                        @inbounds x[idx] = parse(Float64, fields[1])
                        @inbounds y[idx] = parse(Float64, fields[2])
                end
        end

        # Field
        line = readline(fh)
        fields = split(line)
        @assert fields[1] == "POINT_DATA"
        @assert parse(Int64, fields[2]) == numpts
        line = readline(fh)
        @assert line == "FIELD scalar 1"
        line = readline(fh)
        fields = split(line)
        @assert fields[1] == "data"
        @assert fields[2] == "1"
        @assert parse(Int64, fields[3]) == numpts
        @assert fields[4] == "float"

        for j in 1:ny
                for i in 1:nx
                        line = readline(fh)
                        idx = j + (i - 1) * ny
                        @inbounds field[idx] = parse(Float64, line)
                end
        end

        close(fh)
        return x, y, field

end


end
