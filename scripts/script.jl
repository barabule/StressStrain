import Pkg; Pkg.activate("scripts")
#### ]dev ./../StressStrain ###
using StressStrain
SS = StressStrain



# fn = "assets/toein_pp.txt"
# data = SS.read_stress_strain_data(fn;delim='\t',
#         strain_multiplier = 1e-2,
#         )


# fn = "assets/DC04_QS.csv"
# data = SS.read_stress_strain_data(fn; delim = ',',
#                         skipstart = 1,
#                         )

fn= "assets/toein_pp.txt"
data = SS.read_stress_strain_data(fn; delim = '\t',
                            strain_multiplier = 1e-2)

# fn = "assets/noisy.txt"
# data = SS.read_stress_strain_data(fn; delim = '\t',
#                         skipstart = 0,
#                         )

SS.main(data)