@testset "Downscale" begin
  @test !isrevertible(Downscale(2, 2))

  grid = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(10, 10))
  tgrid = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(20, 20))
  gtb = georef((a=rand(Float64, 100), b=rand(Int, 100)), grid)
  ngtb = gtb |> Downscale(2, 2)
  @test domain(ngtb) == tgrid
  @test length(ngtb.a) == nelements(domain(ngtb))
  @test length(ngtb.b) == nelements(domain(ngtb))
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
  @test length(ngtb.a) == nelements(domain(ngtb))
  @test length(ngtb.b) == nelements(domain(ngtb))
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

  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  gtb = georef((a=rand(Float64, 100), b=rand(Int, 100)), sgrid)
  ngtb = gtb |> Downscale(2, 2)
  @test domain(ngtb) == tsgrid
  @test length(ngtb.a) == nelements(domain(ngtb))
  @test length(ngtb.b) == nelements(domain(ngtb))
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

  grid = CartesianGrid(3, 3, 3)
  tgrid = CartesianGrid(minimum(grid), maximum(grid), dims=(6, 6, 6))
  gtb = georef((a=rand(Float64, 27), b=rand(Int, 27)), grid)
  ngtb = gtb |> Downscale(2, 2, 2)
  @test domain(ngtb) == tgrid
  @test length(ngtb.a) == nelements(domain(ngtb))
  @test length(ngtb.b) == nelements(domain(ngtb))
  @test ngtb[(1, 1, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 2, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 1, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 2, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 1, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 2, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 1, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 2, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(5, 5, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 6, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 5, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 6, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 5, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 6, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 5, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 6, 6), :b] == gtb[(3, 3, 3), :b]

  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  gtb = georef((a=rand(Float64, 27), b=rand(Int, 27)), rgrid)
  ngtb = gtb |> Downscale(2, 2, 2)
  @test domain(ngtb) == trgrid
  @test length(ngtb.a) == nelements(domain(ngtb))
  @test length(ngtb.b) == nelements(domain(ngtb))
  @test ngtb[(1, 1, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 2, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 1, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 2, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 1, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 2, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 1, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 2, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(5, 5, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 6, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 5, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 6, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 5, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 6, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 5, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 6, 6), :b] == gtb[(3, 3, 3), :b]

  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  gtb = georef((a=rand(Float64, 27), b=rand(Int, 27)), sgrid)
  ngtb = gtb |> Downscale(2, 2, 2)
  @test domain(ngtb) == tsgrid
  @test length(ngtb.a) == nelements(domain(ngtb))
  @test length(ngtb.b) == nelements(domain(ngtb))
  @test ngtb[(1, 1, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 2, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 1, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 2, 1), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 1, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(1, 2, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 1, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(2, 2, 2), :a] == gtb[(1, 1, 1), :a]
  @test ngtb[(5, 5, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 6, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 5, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 6, 5), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 5, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(5, 6, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 5, 6), :b] == gtb[(3, 3, 3), :b]
  @test ngtb[(6, 6, 6), :b] == gtb[(3, 3, 3), :b]
end
