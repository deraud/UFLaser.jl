using Random, Plots

# Generate a random binary matrix (0s and 1s)
n, m = 20, 20  # Matrix size
Random.seed!(42)  # For reproducibility
matrix = rand(0:1, n, m)

# Define colormap: 0 (vapor) -> blue, 1 (solid) -> red
cmap = [:blue, :red]

# Create X and Y grid
x = 1:m
y = 1:n

# Plot filled contour
contourf(x, y, matrix, 
    fill=true, 
    color=cmap, 
    linewidth=0,  
    colorbar=false,
    title="Phase Distribution")

# Add labels manually
annotate!(m - 2, n - 2, text("Vapor (0)", 10, :blue, :white))
annotate!(m -2, n - 4, text("Solid (1)", 10, :red, :white))
