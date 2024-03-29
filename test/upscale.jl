@testset "Upscale" begin
  @test !isrevertible(Upscale(2, 2))

  grid = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(20, 20))
  tgrid = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(10, 10))
  gtb = georef((a=rand(Float64, 400), b=rand(Int, 400)), grid)
  ngtb = gtb |> Upscale(2, 2)
  @test domain(ngtb) == tgrid
  @test ngtb[(1, 1), :a] == mean(gtb[(1:2, 1:2), :a])
  @test ngtb[(1, 10), :a] == mean(gtb[(1:2, 19:20), :a])
  @test ngtb[(10, 1), :a] == mean(gtb[(19:20, 1:2), :a])
  @test ngtb[(10, 10), :a] == mean(gtb[(19:20, 19:20), :a])
  @test ngtb[(1, 1), :b] == first(gtb[(1:2, 1:2), :b])
  @test ngtb[(1, 10), :b] == first(gtb[(1:2, 19:20), :b])
  @test ngtb[(10, 1), :b] == first(gtb[(19:20, 1:2), :b])
  @test ngtb[(10, 10), :b] == first(gtb[(19:20, 19:20), :b])

  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  gtb = georef((a=rand(Float64, 400), b=rand(Int, 400)), rgrid)
  ngtb = gtb |> Upscale(2, 2)
  @test domain(ngtb) == trgrid
  @test ngtb[(1, 1), :a] == mean(gtb[(1:2, 1:2), :a])
  @test ngtb[(1, 10), :a] == mean(gtb[(1:2, 19:20), :a])
  @test ngtb[(10, 1), :a] == mean(gtb[(19:20, 1:2), :a])
  @test ngtb[(10, 10), :a] == mean(gtb[(19:20, 19:20), :a])
  @test ngtb[(1, 1), :b] == first(gtb[(1:2, 1:2), :b])
  @test ngtb[(1, 10), :b] == first(gtb[(1:2, 19:20), :b])
  @test ngtb[(10, 1), :b] == first(gtb[(19:20, 1:2), :b])
  @test ngtb[(10, 10), :b] == first(gtb[(19:20, 19:20), :b])

  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  gtb = georef((a=rand(Float64, 400), b=rand(Int, 400)), sgrid)
  ngtb = gtb |> Upscale(2, 2)
  @test domain(ngtb) == tsgrid
  @test ngtb[(1, 1), :a] == mean(gtb[(1:2, 1:2), :a])
  @test ngtb[(1, 10), :a] == mean(gtb[(1:2, 19:20), :a])
  @test ngtb[(10, 1), :a] == mean(gtb[(19:20, 1:2), :a])
  @test ngtb[(10, 10), :a] == mean(gtb[(19:20, 19:20), :a])
  @test ngtb[(1, 1), :b] == first(gtb[(1:2, 1:2), :b])
  @test ngtb[(1, 10), :b] == first(gtb[(1:2, 19:20), :b])
  @test ngtb[(10, 1), :b] == first(gtb[(19:20, 1:2), :b])
  @test ngtb[(10, 10), :b] == first(gtb[(19:20, 19:20), :b])

  tgrid = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(10, 5))
  gtb = georef((a=rand(Float64, 400), b=rand(Int, 400)), grid)
  ngtb = gtb |> Upscale(2, 4)
  @test domain(ngtb) == tgrid
  @test ngtb[(1, 1), :a] == mean(gtb[(1:2, 1:4), :a])
  @test ngtb[(1, 5), :a] == mean(gtb[(1:2, 17:20), :a])
  @test ngtb[(10, 1), :a] == mean(gtb[(19:20, 1:4), :a])
  @test ngtb[(10, 5), :a] == mean(gtb[(19:20, 17:20), :a])
  @test ngtb[(1, 1), :b] == first(gtb[(1:2, 1:4), :b])
  @test ngtb[(1, 5), :b] == first(gtb[(1:2, 17:20), :b])
  @test ngtb[(10, 1), :b] == first(gtb[(19:20, 1:4), :b])
  @test ngtb[(10, 5), :b] == first(gtb[(19:20, 17:20), :b])
end
