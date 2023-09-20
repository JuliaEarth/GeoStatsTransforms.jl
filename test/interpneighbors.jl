@testset "InterpolateNeighbors" begin
  @test !isrevertible(InterpolateNeighbors(CartesianGrid(2, 2)))

  pset = PointSet(rand(Point2, 3))
  gtb = georef((a=[1, 2, 3], b=[4, 5, 6]), pset)
  ngtb = gtb |> InterpolateNeighbors(pset, maxneighbors=3)
  @test ngtb.a == gtb.a
  @test ngtb.b == gtb.b
  @test ngtb.geometry == pset

  gtb = georef((; z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
  grid = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
  linds = LinearIndices(size(grid))
  variogram = GaussianVariogram(range=35.0, nugget=0.0)

  Random.seed!(2021)
  ngtb = gtb |> InterpolateNeighbors(grid, :z => Kriging(variogram), maxneighbors=3)
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)

  Random.seed!(2021)
  ngtb = gtb |> InterpolateNeighbors(grid, :z => Kriging(variogram), maxneighbors=3)
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)

  Random.seed!(2021)
  ngtb = gtb |> InterpolateNeighbors(grid, :z => Kriging(variogram), maxneighbors=3, neighborhood=MetricBall(100.0))
  @test isapprox(ngtb.z[linds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(ngtb.z[linds[75, 50]], 1.0, atol=1e-3)
end
