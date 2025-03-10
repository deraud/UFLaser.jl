using Interpolations

# Given data points

#n
T = [300, 15000, 30000]
n = [3.58, 2.94, 2.55]

#k
T = [300, 15000, 30000]
k = [4.62, 3.31, 2.91]

# Create a linear interpolation function
interp = LinearInterpolation(T, n, extrapolation_bc=Flat())

# Function to interpolate or extrapolate with constant bounds
function interpolate_value(x_val)
    if x_val < 300
        return n[1]
    elseif x_val > 30000
        return n[end]
    else
        return interp(x_val)
    end
end

# Example usage

@time for i in 1:1250000
    interpolate_value(5000)  # Interpolation example
end

