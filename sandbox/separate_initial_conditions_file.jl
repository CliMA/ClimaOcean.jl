using JLD2
using Printf

filename = "initial_conditions_T_S_360x150x48x12.jld2"

file = jldopen(filename)
T = file["T"]
S = file["S"]
close(file)

for month = 1:12
    Tm = T[:, :, :, month]
    Sm = S[:, :, :, month]
    savename = @sprintf("initial_conditions_month_%02d_360_150_48.jld2",
                        month)

    @save savename T=Tm S=Sm
end

