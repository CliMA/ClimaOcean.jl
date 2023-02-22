using JLD2

#####
##### One degree boundary conditions
#####

onedegreepath = "surface_boundary_conditions_12_months_360_150.jld2"

file = jldopen(onedegreepath)
east_momentum_flux  = - file["τˣ"]
north_momentum_flux = - file["τʸ"]
heat_flux           = - file["Qᶠ"]
salinity_flux       = - file["Sᶠ"]
surface_temperature = + file["Tₛ"]
surface_salinity    = + file["Sₛ"]
close(file)

new_onedegreepath = "near_global_boundary_conditions_360_150.jld2"
isfile(new_onedegreepath) && rm(new_onedegreepath)
file = jldopen(new_onedegreepath, "a+")
file["east_momentum_flux"]  = east_momentum_flux  
file["north_momentum_flux"] = north_momentum_flux 
file["heat_flux"]           = heat_flux           
file["salinity_flux"]       = salinity_flux       
file["surface_temperature"] = surface_temperature 
file["surface_salinity"]    = surface_salinity    
close(file)

#####
##### Quarter degree boundary conditions
#####

temppath = "temp-1440x600-latitude-75.jld2"
saltpath = "salt-1440x600-latitude-75.jld2"
tauxpath = "tau_x-1440x600-latitude-75.jld2"
tauypath = "tau_y-1440x600-latitude-75.jld2"

file = jldopen(temppath)
surface_temperature = file["field"]
close(file)

file = jldopen(saltpath)
surface_salinity = file["field"]
close(file)

file = jldopen(tauxpath)
east_momentum_flux = - file["field"]
close(file)

file = jldopen(tauypath)
north_momentum_flux = - file["field"]
close(file)

new_quarterdegreepath = "near_global_boundary_conditions_1440_600.jld2"
isfile(new_quarterdegreepath) && rm(new_quarterdegreepath)
file = jldopen(new_quarterdegreepath, "a+")
file["east_momentum_flux"]  = east_momentum_flux  
file["north_momentum_flux"] = north_momentum_flux 
file["surface_temperature"] = surface_temperature 
file["surface_salinity"]    = surface_salinity    
close(file)

path = "near_global_east_momentum_flux_1440_600.jld2"
isfile(path) && rm(path)
file = jldopen(path, "a+")
file["east_momentum_flux"] = east_momentum_flux  
close(file)

path = "near_global_north_momentum_flux_1440_600.jld2"
isfile(path) && rm(path)
file = jldopen(path, "a+")
file["north_momentum_flux"] = north_momentum_flux  
close(file)

path = "near_global_surface_temperature_1440_600.jld2"
isfile(path) && rm(path)
file = jldopen(path, "a+")
file["surface_temperature"] = surface_temperature  
close(file)

path = "near_global_surface_salinity_1440_600.jld2"
isfile(path) && rm(path)
file = jldopen(path, "a+")
file["surface_salinity"] = surface_salinity  
close(file)

