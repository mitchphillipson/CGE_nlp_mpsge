using MPSGE


begin
    MP = MPSGEModel()

    @sectors(MP, begin
        X
        Q
    end)

    @commodities(MP, begin
        PFX 
        PX
        PQ
        PDT
    end)

    @consumers(MP, begin
        world
        Y
    end)


    @production(MP, X, [s=.2,t=0], begin
        @output(PX, 100, t)
        @input(PFX, 25, s)
        @input(PDT, 75, s)
    end)


    @production(MP, Q, [s=0,t=.2], begin
        @output(PFX, 25, t)
        @output(PDT, 75, t)
        @input(PQ, 100, s)
    end)


    @demand(MP, world, begin
        @final_demand(PFX, 25) #export
        @endowment(PFX, 25) #import
    end)

    @demand(MP, Y, begin
        @final_demand(PFX, 0)
        @final_demand(PX, 100)
        @endowment(PQ, 100)
    end)
end

fix(MP[:PQ], 1)

solve!(MP)


value(MP[:PFX])