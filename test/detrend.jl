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

  rng = MersenneTwister(42)
  # constant trend
  d = georef((z=rand(rng, 100),), CartesianGrid(100))
  z̄ = GeoStatsTransforms.trend(d, :z).z
  @test all(abs.(diff(z̄)) .< 0.01)

  # linear trend
  μ = range(0, stop=1, length=100)
  ϵ = 0.1rand(rng, 100)
  d = georef((z=μ + ϵ,), CartesianGrid(100))
  z̄ = GeoStatsTransforms.trend(d, :z).z
  @test all([abs(z̄[i] - μ[i]) < 0.1 for i in 1:length(z̄)])

  # quadratic trend
  r = range(-1, stop=1, length=100)
  μ = [x^2 + y^2 for x in r, y in r]
  ϵ = 0.1rand(rng, 100, 100)
  d = georef((z=μ + ϵ,))
  d̄ = GeoStatsTransforms.trend(d, :z, degree=2)
  z̄ = reshape(d̄.z, 100, 100)
  @test all([abs(z̄[i] - μ[i]) < 0.1 for i in 1:length(z̄)])

  d = georef((x=rand(rng, 10), y=rand(rng, 10)), rand(rng, Point, 10))
  𝒯 = d |> GeoStatsTransforms.trend |> values
  s = Tables.schema(𝒯)
  @test s.names == (:x, :y)
  @test s.types == (Float64, Float64)
end
