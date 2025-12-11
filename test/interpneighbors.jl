@testset "InterpolateNeighbors" begin
  @test !isrevertible(InterpolateNeighbors(CartesianGrid(2, 2)))

  pts = rand(Point, 3)
  gtb = georef((; z=[1, 2, 3]), pts)
  ngtb = gtb |> InterpolateNeighbors(pts; model=IDW(), maxneighbors=3)
  @test ngtb.z == gtb.z
  @test ngtb.geometry == gtb.geometry

  gtb = georef((; z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
  grid = CartesianGrid((0.5, 0.5), (100.5, 100.5), dims=(100, 100))
  inds = LinearIndices(size(grid))
  model = Kriging(GaussianVariogram(range=35.0))
  ngtb = gtb |> InterpolateNeighbors(grid; model=model)
  @test isapprox(ngtb.z[inds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[inds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[inds[75, 50]], 1.0, atol=1e-3)
  ngtb = gtb |> InterpolateNeighbors(grid; model=model, neighborhood=MetricBall(100.0))
  @test isapprox(ngtb.z[inds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[inds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[inds[75, 50]], 1.0, atol=1e-3)

  # units
  gtb = georef((; T=[1.0, 0.0, 1.0] * u"K"), rand(Point, 3))
  grid = CartesianGrid(5, 5, 5)
  ngtb = gtb |> InterpolateNeighbors(grid; model=IDW())
  @test unit(eltype(ngtb.T)) == u"K"

  # affine units
  gtb = georef((; T=[-272.15, -273.15, -272.15] * u"°C"), rand(Point, 3))
  ngtb = gtb |> InterpolateNeighbors(grid; model=IDW())
  @test unit(eltype(ngtb.T)) == u"K"

  # default model is NN
  pts = rand(Point, 3)
  gtb = georef((; z=["a", "b", "c"]), pts)
  ngtb = gtb |> InterpolateNeighbors(pts)
  @test ngtb == gtb
end
