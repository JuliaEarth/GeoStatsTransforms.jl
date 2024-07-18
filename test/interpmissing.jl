@testset "InterpolateMissing" begin
  @test !isrevertible(InterpolateMissing())

  grid = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
  linds = LinearIndices(size(grid))
  z = Vector{Union{Missing,Float64}}(missing, 10000)
  z[linds[25, 25]] = 1.0
  z[linds[50, 75]] = 0.0
  z[linds[75, 50]] = 1.0
  gtb = georef((; z), grid)
  variogram = GaussianVariogram(range=35.0, nugget=0.0)

  Random.seed!(2021)
  ngtb = gtb |> InterpolateMissing(:z => Kriging(variogram), maxneighbors=3)
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)

  Random.seed!(2021)
  ngtb = gtb |> InterpolateMissing(:z => Kriging(variogram), maxneighbors=3)
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)

  Random.seed!(2021)
  ngtb = gtb |> InterpolateMissing(:z => Kriging(variogram), maxneighbors=3, neighborhood=MetricBall(100.0))
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)

  # units
  grid = CartesianGrid(5, 5)
  gtb = georef((; T=shuffle([[1.0, 0.0, 1.0]; fill(missing, 25 - 3)]) * u"K"), grid)
  ngtb = gtb |> InterpolateMissing(IDW())
  @test unit(eltype(ngtb.T)) == u"K"

  # affine units
  gtb = georef((; T=shuffle([[-272.15, -273.15, -272.15]; fill(missing, 25 - 3)]) * u"°C"), grid)
  ngtb = gtb |> InterpolateMissing(IDW())
  @test unit(eltype(ngtb.T)) == u"K"

  # default model is NN
  Random.seed!(2021)
  pset = PointSet(rand(Point, 5))
  gtb = georef((; z=[1.0, missing, 2.0, missing, 3.0], c=[missing, "a", "b", "c", missing]), pset)
  ngtb = gtb |> InterpolateMissing()
  @test ngtb.z ⊆ [1.0, 2.0, 3.0]
  @test ngtb.c ⊆ ["a", "b", "c"]
end
