@testset "Transfer" begin
  @test !isrevertible(Transfer(CartesianGrid(10, 10)))

  pts1 = Point.([(1, 1), (7, 1), (4, 4)])
  pts2 = Point.([(5, 4), (3, 4), (0, 1), (7, 0), (7, 2)])
  gtb = georef((a=rand(Float64, 3), b=rand(Int, 3)), pts1)
  ngtb = gtb |> Transfer(pts2)
  @test domain(ngtb) == PointSet(pts2)
  @test length(ngtb.a) == nelements(domain(ngtb))
  @test length(ngtb.b) == nelements(domain(ngtb))
  @test ngtb.a[1] == gtb.a[3]
  @test ngtb.a[2] == gtb.a[3]
  @test ngtb.a[3] == gtb.a[1]
  @test ngtb.a[4] == gtb.a[2]
  @test ngtb.a[5] == gtb.a[2]
  @test ngtb.b[1] == gtb.b[3]
  @test ngtb.b[2] == gtb.b[3]
  @test ngtb.b[3] == gtb.b[1]
  @test ngtb.b[4] == gtb.b[2]
  @test ngtb.b[5] == gtb.b[2]

  grid1 = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(10, 10))
  grid2 = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(20, 20))
  gtb = georef((a=rand(Float64, 100), b=rand(Int, 100)), grid1)
  ngtb = gtb |> Transfer(grid2)
  @test domain(ngtb) == grid2
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
end
