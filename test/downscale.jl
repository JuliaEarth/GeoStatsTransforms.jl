@testset "Downscale" begin
  @test !isrevertible(Downscale(2, 2))
  
  grid = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(10, 10))
  tgrid = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(20, 20))
  gtb = georef((a=rand(Float64, 100), b=rand(Int, 100)), grid)
  ngtb = gtb |> Downscale(2, 2)
  @test domain(ngtb) == tgrid
  @test ngtb[(1, 1), :a] == gtb[(1, 1), :a]
  @test ngtb[(1, 2), :a] == gtb[(1, 1), :a]
  @test ngtb[(2, 1), :a] == gtb[(1, 1), :a]
  @test ngtb[(2, 2), :a] == gtb[(1, 1), :a]
  @test ngtb[(1, 19), :a] == gtb[(1, 10), :a]
  @test ngtb[(1, 20), :a] == gtb[(1, 10), :a]
  @test ngtb[(2, 19), :a] == gtb[(1, 10), :a]
  @test ngtb[(2, 20), :a] == gtb[(1, 10), :a]
  @test ngtb[(19, 1), :b] == gtb[(10, 1), :b]
  @test ngtb[(19, 2), :b] == gtb[(10, 1), :b]
  @test ngtb[(20, 1), :b] == gtb[(10, 1), :b]
  @test ngtb[(20, 2), :b] == gtb[(10, 1), :b]
  @test ngtb[(19, 19), :b] == gtb[(10, 10), :b]
  @test ngtb[(19, 20), :b] == gtb[(10, 10), :b]
  @test ngtb[(20, 19), :b] == gtb[(10, 10), :b]
  @test ngtb[(20, 20), :b] == gtb[(10, 10), :b]

  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  gtb = georef((a=rand(Float64, 100), b=rand(Int, 100)), rgrid)
  ngtb = gtb |> Downscale(2, 2)
  @test domain(ngtb) == trgrid
  @test ngtb[(1, 1), :a] == gtb[(1, 1), :a]
  @test ngtb[(1, 2), :a] == gtb[(1, 1), :a]
  @test ngtb[(2, 1), :a] == gtb[(1, 1), :a]
  @test ngtb[(2, 2), :a] == gtb[(1, 1), :a]
  @test ngtb[(1, 19), :a] == gtb[(1, 10), :a]
  @test ngtb[(1, 20), :a] == gtb[(1, 10), :a]
  @test ngtb[(2, 19), :a] == gtb[(1, 10), :a]
  @test ngtb[(2, 20), :a] == gtb[(1, 10), :a]
  @test ngtb[(19, 1), :b] == gtb[(10, 1), :b]
  @test ngtb[(19, 2), :b] == gtb[(10, 1), :b]
  @test ngtb[(20, 1), :b] == gtb[(10, 1), :b]
  @test ngtb[(20, 2), :b] == gtb[(10, 1), :b]
  @test ngtb[(19, 19), :b] == gtb[(10, 10), :b]
  @test ngtb[(19, 20), :b] == gtb[(10, 10), :b]
  @test ngtb[(20, 19), :b] == gtb[(10, 10), :b]
  @test ngtb[(20, 20), :b] == gtb[(10, 10), :b]
end
