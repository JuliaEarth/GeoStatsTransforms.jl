@testset "Detrend" begin
  rng = MersenneTwister(42)
  l = range(-1, stop=1, length=100)
  μ = [x^2 + y^2 for x in l, y in l]
  ϵ = 0.1rand(rng, 100, 100)
  d = georef((z=μ + ϵ, w=rand(100, 100)))
  p = Detrend(:z, degree=2)
  n, c = apply(p, d)
  r = revert(p, n, c)
  D = Tables.matrix(values(d))
  R = Tables.matrix(values(r))
  @test isapprox(D, R, atol=1e-6)
end
