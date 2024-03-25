@testset "Aggregate" begin
  @test !isrevertible(Aggregate(CartesianGrid(10, 10)))

  pts1 = Point2[(5, 4), (3, 4), (0, 1), (7, 0), (7, 2)]
  pts2 = Point2[(1, 1), (7, 1), (4, 4)]
  gtb = georef((a=rand(Float64, 5), b=rand(Int, 5)), pts1)
  ngtb = gtb |> Aggregate(pts2)
  @test domain(ngtb) == PointSet(pts2)
  @test ngtb.a[1] == gtb.a[3]
  @test ngtb.a[2] == mean(gtb.a[[4, 5]])
  @test ngtb.a[3] == mean(gtb.a[[1, 2]])
  @test ngtb.b[1] == gtb.b[3]
  @test ngtb.b[2] == first(gtb.b[[4, 5]])
  @test ngtb.b[3] == first(gtb.b[[1, 2]])

  ngtb = gtb |> Aggregate(pts2, :a => median, :b => last)
  @test domain(ngtb) == PointSet(pts2)
  @test ngtb.a[1] == gtb.a[3]
  @test ngtb.a[2] == median(gtb.a[[4, 5]])
  @test ngtb.a[3] == median(gtb.a[[1, 2]])
  @test ngtb.b[1] == gtb.b[3]
  @test ngtb.b[2] == last(gtb.b[[4, 5]])
  @test ngtb.b[3] == last(gtb.b[[1, 2]])

  grid1 = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(20, 20))
  grid2 = CartesianGrid((0.0, 0.0), (10.0, 10.0), dims=(10, 10))
  gtb = georef((a=rand(Float64, 400), b=rand(Int, 400)), grid1)
  ngtb = gtb |> Aggregate(grid2)
  @test domain(ngtb) == grid2
  @test ngtb[(1, 1), :a] == mean(gtb[(1:2, 1:2), :a])
  @test ngtb[(1, 10), :a] == mean(gtb[(1:2, 19:20), :a])
  @test ngtb[(10, 1), :a] == mean(gtb[(19:20, 1:2), :a])
  @test ngtb[(10, 10), :a] == mean(gtb[(19:20, 19:20), :a])
  @test ngtb[(1, 1), :b] == first(gtb[(1:2, 1:2), :b])
  @test ngtb[(1, 10), :b] == first(gtb[(1:2, 19:20), :b])
  @test ngtb[(10, 1), :b] == first(gtb[(19:20, 1:2), :b])
  @test ngtb[(10, 10), :b] == first(gtb[(19:20, 19:20), :b])

  ngtb = gtb |> Aggregate(grid2, :a => median, :b => last)
  @test domain(ngtb) == grid2
  @test ngtb[(1, 1), :a] == median(gtb[(1:2, 1:2), :a])
  @test ngtb[(1, 10), :a] == median(gtb[(1:2, 19:20), :a])
  @test ngtb[(10, 1), :a] == median(gtb[(19:20, 1:2), :a])
  @test ngtb[(10, 10), :a] == median(gtb[(19:20, 19:20), :a])
  @test ngtb[(1, 1), :b] == last(gtb[(1:2, 1:2), :b])
  @test ngtb[(1, 10), :b] == last(gtb[(1:2, 19:20), :b])
  @test ngtb[(10, 1), :b] == last(gtb[(19:20, 1:2), :b])
  @test ngtb[(10, 10), :b] == last(gtb[(19:20, 19:20), :b])
end
