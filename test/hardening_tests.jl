
@testset "Data Prep" verbose= true begin
    E = round(1e5 * rand() + 1e5; sigdigits = 2)

    strain = LinRange(0, 1e-2, 10)
    stress = E .* strain

    data = (;strain , stress)

    @test SS.get_modulus(data; sigdigits = 2) ≈ E

    Ebrk = SS.bracket_modulus(data)
    @test Ebrk[1] ≈ Ebrk[2] ≈ E


    circ = SS.true_to_engineering(SS.engineering_to_true(data))
    @test all(circ.strain .≈ data.strain)
    @test all(circ.stress .≈ data.stress)

    circ = SS.engineering_to_true(SS.true_to_engineering(data))
    @test all(circ.strain .≈ data.strain)
    @test all(circ.stress .≈ data.stress)


    Etan = 1e-2 * E
    sy = 200
    ey = sy/ E
    data = (;strain = [0.0, 0.001, ey, 0.1, 0.2],
            stress = [0.0, E * 1e-3, sy, Etan * (0.1 - ey), Etan * (0.2 - ey)])

    hardening = SS.get_hardening_portion(data, E)
    @test length(hardening.strain) == 2

    Ebrk2 = SS.bracket_modulus(data;sigdigits =2)

    @test Ebrk2.Emax ≈ E
    @test Ebrk2.Emin <= E <= Ebrk2.Emax


    offset = 2e-3

    hardening = SS.get_hardening_portion(data, E; offset)
    @test issorted(hardening.strain)
    @test begin
        val = true
        etotal = view(data.strain, (length(data.strain)-length(hardening.strain)+1) : length(data.strain))
        for (epl, sef, et)  in zip(hardening.strain, hardening.stress, etotal)
            eel = sef / E
            if !(et == (eel + epl))
                val = false
                @info " ϵt = $(et), ϵpl = $epl, ϵel = $eel"
                break
            end
        end 
        val
    end
    # @test (last(hardening.strain)  + (last(hardening.stress) / E) - last(data.strain)) ≈ 0
    @test hardening.stress[1] == data.stress[4]


    #toein stuff
    cut = 1e-2
    toein = SS.toein_compensate(data; cut)
    @test issorted(toein.strain)
    @test first(toein.strain) ≈ 0 && first(toein.stress) ≈ 0 #after toein operation, stress strain should start at 0

    nl = length(toein.stress)
    @test all(data.stress[end-nl+2:end] .== toein.stress[2:end])  #should not modify stress values after the cut
    @test all(data.strain[end-nl+2:end] .>= cut) # strain after cut
    ncut = length(data.strain) - nl
    @test all(data.strain[1:ncut] .< cut) #strain before cut should be < cut

end
        