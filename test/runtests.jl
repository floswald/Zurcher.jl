using Test, Zurcher


@testset "transition matrix" begin
    θ = [0.1,0.2,0.6,0.1]
    n = 37
    p = Zurcher.make_trans(θ,n)
    @test size(p) == (n,n)
    @test all(sum(p,dims = 2) .== 1.0)
    @test p[2,2:(2+length(θ)-1)] == θ
end

@testset "Building Harold" begin
    n = 105
    h = Zurcher.Harold(n = n)
    @test h.n == n
    @test length(h.mileage) == n

    ccp = Zurcher.ccp(h,ones(n))
    @test all((ccp .>= 0) .| (ccp .<= 1.0))

end