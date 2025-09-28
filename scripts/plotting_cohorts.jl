
using CairoMakie, ColorSchemes
using Connectomes
using ComponentArrays
using CSV
using DataFrames
using GLMakie; GLMakie.activate!()

# Load the data from the updated CSV - here is example for homogeneous clearance simulation outputs 
df = CSV.read("uni_output.csv", DataFrame)

# Load the connectome and filter for the right hemisphere cortex
_c = Connectome(Connectomes.connectome_path())
cortex = filter(x -> get_lobe(x) !== "subcortex" && get_hemisphere(x) == "right", _c.parc)
c = filter(slice(_c, cortex))
nodes = get_node_id.(cortex)  # Get nodes corresponding to the right hemisphere connectome

# Set the color range
simMinval = 0
simMaxval = 1

cmap=ColorSchemes.viridis

f_cbar = Figure(size=(50, 100))
 Colorbar(f_cbar[1, 1], 
 limits=(simMinval, simMaxval), 
 colormap=cmap, 
 label="tp", 
 labelsize=10,    
 ticklabelsize=6) 
 save("plots/colorbar.png", f_cbar, px_per_unit=20)

# Function to create and save a plot for a given output (front and back views)
function create_patient_plot(df, timepoint)
    output_values = zeros(length(nodes), 1)
    
    # Set initial conditions 
    for (i, row) in enumerate(eachrow(df))
        value = row[timepoint]
        output_values[i, 1] = value
    end

    print(output_values)

    # Cap the values at the specified bounds
    capped_values = clamp.(output_values, simMinval, simMaxval)

    ics = ComponentArray(v=capped_values)

    # Normalize the data for plotting
    max_norm(x) = x ./ maximum(x)

    # Front View
    f_front = Figure(size=(300, 220))  
    ax = Axis3(f_front[1,1], 
    aspect = :data, 
    azimuth = 1.0pi, 
    elevation=0.0pi,
    protrusions=(1.0,1.0,1.0,1.0))
    hidedecorations!(ax)
    hidespines!(ax)
    plot_roi!(nodes, max_norm(ics.v[:, 1]), cmap)
    save("plots/uni_ave_output_$(timepoint)_front.png", f_front, px_per_unit=20)

    # Back View
    f_back = Figure(size=(300, 220))  
    ax = Axis3(f_back[1,1], 
    aspect = :data, 
    azimuth = 0.0pi, 
    elevation=0.0pi,
    protrusions=(1.0,1.0,1.0,1.0))
    hidedecorations!(ax)
    hidespines!(ax)
    plot_roi!(nodes, max_norm(ics.v[:, 1]), cmap)
    save("plots/uni_ave_output_$(timepoint)_back.png", f_back, px_per_unit=20)
end

# Loop over each time point in the DataFrame and create front and back plots for each
for timepoint in names(df)[3:end]  
    create_patient_plot(df, timepoint)
end
