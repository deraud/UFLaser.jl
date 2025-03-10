include("./VOF.jl")

module TTM_VOF_Coupling

using ..FreeSurfaceLBM
using ..Structs
using ..Constants
using ..Thermal



export initializeVOFFromTTM
export updateVOFFromTTM
export mapCellTypes

# Temperature thresholds for phase changes
const MELTING_TEMP = 1337.0  # Melting point of gold
const VAPORIZATION_TEMP = 3129.0  # Boiling point of gold

function initializeVOFFromTTM(ttm_lattice::Lattice, vof_domain::LBMDomain)
    # Initialize VOF domain based on initial TTM state
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            # Map TTM grid to VOF grid (may need interpolation if grids differ)
            ttm_i = min(max(round(Int, ix * Nx / vof_domain.nx), 1), Nx)
            ttm_j = min(max(round(Int, iy * Ny / vof_domain.ny), 1), Ny)
            
            # Set cell type based on temperature
            vof_domain.cell_type[ix, iy] = mapCellType(ttm_lattice.Tₗ[ttm_i, ttm_j])
            
            # Set initial density and velocity
            if vof_domain.cell_type[ix, iy] & CellType.FLUID != 0
                vof_domain.rho[ix, iy] = 1.0
                vof_domain.mass[ix, iy] = 1.0
            elseif vof_domain.cell_type[ix, iy] & CellType.INTERFACE != 0
                vof_domain.rho[ix, iy] = 1.0
                vof_domain.mass[ix, iy] = 0.5
            else
                vof_domain.rho[ix, iy] = 0.0
                vof_domain.mass[ix, iy] = 0.0
            end
            
            # Initialize distribution functions with equilibrium
            FreeSurfaceLBM.set_equi!(vof_domain, ix, iy, vof_domain.rho[ix, iy], vof_domain.u[:, ix, iy])
        end
    end
    
    # Copy equilibrium to distribution functions
    vof_domain.fin .= vof_domain.equi
    
    return vof_domain
end

function mapCellType(temperature::Float64)
    if temperature >= VAPORIZATION_TEMP
        return CellType.GAS
    elseif temperature >= MELTING_TEMP
        return CellType.FLUID
    else
        return CellType.OBSTACLE  # Solid material
    end
end

function updateVOFFromTTM(ttm_lattice::Lattice, vof_domain::LBMDomain)
    # Update VOF domain based on TTM state
    
    # 1. Apply thermal forces from temperature gradients
    applyThermalForces!(ttm_lattice, vof_domain)
    
    # 2. Update cell types based on temperature
    updateCellTypesFromTemp!(ttm_lattice, vof_domain)
    
    # 3. Map ablation front from TTM to VOF
    mapAblationFront!(ttm_lattice, vof_domain)
    
    return vof_domain
end

function applyThermalForces!(ttm_lattice::Lattice, vof_domain::LBMDomain)
    # Thermal expansion coefficient (example value, should be material-specific)
    thermal_expansion_coef = 1.4e-5
    
    # Calculate thermal forces for each cell
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            # Skip cells that aren't fluid or interface
            if (vof_domain.cell_type[ix, iy] & (CellType.FLUID | CellType.INTERFACE)) == 0
                continue
            end
            
            # Map VOF coordinates to TTM coordinates
            ttm_i = min(max(round(Int, ix * Nx / vof_domain.nx), 1), Nx)
            ttm_j = min(max(round(Int, iy * Ny / vof_domain.ny), 1), Ny)
            
            # Calculate temperature gradients (central difference)
            if ttm_i > 1 && ttm_i < Nx
                grad_temp_x = (ttm_lattice.Tₗ[ttm_i+1, ttm_j] - ttm_lattice.Tₗ[ttm_i-1, ttm_j]) / (2*dx)
            else
                grad_temp_x = 0.0
            end
            
            if ttm_j > 1 && ttm_j < Ny
                grad_temp_y = (ttm_lattice.Tₗ[ttm_i, ttm_j+1] - ttm_lattice.Tₗ[ttm_i, ttm_j-1]) / (2*dy)
            else
                grad_temp_y = 0.0
            end
            
            # Convert temperature gradients to forces
            thermal_force_x = -thermal_expansion_coef * grad_temp_x
            thermal_force_y = -thermal_expansion_coef * grad_temp_y
            
            # Add thermal forces to gravity
            vof_domain.gravity[1] = thermal_force_x
            vof_domain.gravity[2] = -0.1 + thermal_force_y  # Combine with existing gravity
        end
    end
end

function updateCellTypesFromTemp!(ttm_lattice::Lattice, vof_domain::LBMDomain)
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            # Map VOF coordinates to TTM coordinates
            ttm_i = min(max(round(Int, ix * Nx / vof_domain.nx), 1), Nx)
            ttm_j = min(max(round(Int, iy * Ny / vof_domain.ny), 1), Ny)
            
            # Skip if the cell is already marked for change
            if (vof_domain.cell_type[ix, iy] & (CellType.TO_FLUID | CellType.TO_GAS)) != 0
                continue
            end
            
            # Update cell types based on temperature
            temperature = ttm_lattice.Tₗ[ttm_i, ttm_j]
            
            if temperature >= VAPORIZATION_TEMP && 
               (vof_domain.cell_type[ix, iy] & (CellType.FLUID | CellType.INTERFACE)) != 0
                vof_domain.cell_type[ix, iy] |= CellType.TO_GAS
            elseif temperature >= MELTING_TEMP && 
                  (vof_domain.cell_type[ix, iy] & CellType.OBSTACLE) != 0
                vof_domain.cell_type[ix, iy] = CellType.INTERFACE
                # Initialize new interface cell
                average_surround!(vof_domain, ix, iy)
            end
        end
    end
end

function mapAblationFront!(ttm_lattice::Lattice, vof_domain::LBMDomain)
    # Use the ablation front calculated in TTM to update VOF cell types
    for j in 1:min(Ny, vof_domain.ny)
        # Map TTM ablation edge to VOF grid
        edge_ix = min(max(round(Int, ttm_lattice.AblationEdge[j] * vof_domain.nx / Nx), 1), vof_domain.nx)
        
        # Mark cells beyond ablation edge as gas
        for ix in 1:edge_ix
            if (vof_domain.cell_type[ix, j] & (CellType.OBSTACLE | CellType.FLUID)) != 0
                vof_domain.cell_type[ix, j] = CellType.INTERFACE
                vof_domain.cell_type[ix, j] |= CellType.TO_GAS
            end
        end
        
        # Mark interface cells at the ablation front
        if edge_ix < vof_domain.nx
            vof_domain.cell_type[edge_ix+1, j] = CellType.INTERFACE
        end
    end
end

# Update relaxation parameter based on temperature
function updateRelaxationParameter!(ttm_lattice::Lattice, vof_domain::LBMDomain)
    for ix in 1:vof_domain.nx
        for iy in 1:vof_domain.ny
            # Map VOF coordinates to TTM coordinates
            ttm_i = min(max(round(Int, ix * Nx / vof_domain.nx), 1), Nx)
            ttm_j = min(max(round(Int, iy * Ny / vof_domain.ny), 1), Ny)
            
            # Get temperature from TTM
            temperature = ttm_lattice.Tₗ[ttm_i, ttm_j]
            
            # Temperature-dependent viscosity (simplified model)
            # For gold: viscosity decreases with increasing temperature
            if temperature > MELTING_TEMP
                # Linear model: viscosity decreases as temperature increases
                viscosity_factor = max(0.1, 1.0 - 0.0005 * (temperature - MELTING_TEMP))
                # Update relaxation parameter (related to viscosity)
                vof_domain.omega = 1.0 / (3.0 * viscosity_factor + 0.5)
            end
        end
    end
end

end # module TTM_VOF_Coupling