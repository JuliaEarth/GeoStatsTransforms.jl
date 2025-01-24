@testset "Interpolate" begin
  @test !isrevertible(Interpolate(CartesianGrid(2, 2)))

  pts = rand(Point, 3)
  gtb = georef((a=[1, 2, 3], b=[4, 5, 6]), pts)
  ngtb = gtb |> Interpolate(pts; model=IDW())
  @test ngtb.a == gtb.a
  @test ngtb.b == gtb.b
  @test ngtb.geometry == gtb.geometry

  gtb = georef((; z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
  grid = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
  inds = LinearIndices(size(grid))
  ngtb = gtb |> Interpolate(grid; model=Kriging(GaussianVariogram(range=35.0)))
  @test isapprox(ngtb.z[inds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[inds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[inds[75, 50]], 1.0, atol=1e-3)

  # units
  gtb = georef((; T=[1.0, 0.0, 1.0] * u"K"), rand(Point, 3))
  grid = CartesianGrid(5, 5, 5)
  ngtb = gtb |> Interpolate(grid)
  @test unit(eltype(ngtb.T)) == u"K"

  # affine units
  gtb = georef((; T=[-272.15, -273.15, -272.15] * u"°C"), rand(Point, 3))
  grid = CartesianGrid(5, 5, 5)
  ngtb = gtb |> Interpolate(grid)
  @test unit(eltype(ngtb.T)) == u"K"

  # default model is NN
  pts = rand(Point, 3)
  gtb = georef((; z=[1.0, 2.0, 3.0], c=["a", "b", "c"]), pts)
  ngtb = gtb |> Interpolate(pts)
  @test ngtb == gtb
end
