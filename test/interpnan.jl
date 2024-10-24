@testset "InterpolateNaN" begin
  @test !isrevertible(InterpolateNaN())

  grid = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
  linds = LinearIndices(size(grid))
  z = fill(NaN, 10000)
  z[linds[25, 25]] = 1.0
  z[linds[50, 75]] = 0.0
  z[linds[75, 50]] = 1.0
  gtb = georef((; z), grid)
  variogram = GaussianVariogram(range=35.0, nugget=0.0)

  ngtb = gtb |> InterpolateNaN(:z => Kriging(variogram), maxneighbors=3)
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)

  ngtb = gtb |> InterpolateNaN(:z => Kriging(variogram), maxneighbors=3)
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)

  ngtb = gtb |> InterpolateNaN(:z => Kriging(variogram), maxneighbors=3, neighborhood=MetricBall(100.0))
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)

  # Kriging with non-unique geometries
  gtb = georef(
    (; z=[1.0, 1.0, NaN, 0.0, 0.0, NaN, 1.0, 1.0]),
    [(1.0, 1.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (3.0, 3.0), (4.0, 4.0), (5.0, 5.0), (5.0, 5.0)]
  )
  variogram = GaussianVariogram(range=35.0, nugget=0.0)
  ngtb = gtb |> InterpolateNaN(:z => Kriging(variogram), maxneighbors=3)
  @test isapprox(ngtb.z[2], 0.25, atol=1e-2)
  @test isapprox(ngtb.z[4], 0.25, atol=1e-2)

  # units
  grid = CartesianGrid(3)
  gtb = georef((; T=[1.0, NaN, 1.0] * u"K"), grid)
  ngtb = gtb |> InterpolateNaN(IDW())
  @test unit(eltype(ngtb.T)) == u"K"

  # affine units
  gtb = georef((; T=[1.0, NaN, 1.0] * u"°C"), grid)
  ngtb = gtb |> InterpolateNaN(IDW())
  @test unit(eltype(ngtb.T)) == u"K"

  # default model is NN
  pset = PointSet(rand(Point, 5))
  gtb = georef((; z=[1.0, NaN, 2.0, NaN, 3.0]), pset)
  ngtb = gtb |> InterpolateNaN()
  @test ngtb.z ⊆ [1.0, 2.0, 3.0]
end
