using CSV, DataFrames

function update_csv!(filename::String, A::Number, B::Number)
    # Check if the file exists
    if isfile(filename)
        df = CSV.read(filename, DataFrame)
    else
        df = DataFrame(A=Int[], B=Float64[])
    end
    
    # Check if A already exists
    if A in df.A
        df[df.A .== A, :B] .= B  # Update existing B value
    else
        push!(df, (A, B))  # Add new row
    end
    
    # Sort by A
    sort!(df, :A)
    
    # Write back to CSV
    CSV.write(filename, df)
end

# Example usage
update_csv!("1D/CSV/dataAblation800ps.csv", trunc(Int,Lx*1e9), J/4.66)