@testset "Aggregate" begin
  @test !isrevertible(Aggregate(CartesianGrid(10, 10)))

  pts1 = Point2[(1, 1), (7, 1), (4, 4)]
  pts2 = Point2[(5, 4), (3, 4), (0, 1), (7, 0), (7, 2)]
  gtb = georef((a=rand(Float64, 5), b=rand(Int, 5)), pts2)
  ngtb = gtb |> Aggregate(pts1)
  @test domain(ngtb) == PointSet(pts1)
  @test ngtb.a[1] == gtb.a[3]
  @test ngtb.a[2] == mean(gtb.a[[4, 5]])
  @test ngtb.a[3] == mean(gtb.a[[1, 2]])
  @test ngtb.b[1] == gtb.b[3]
  @test ngtb.b[2] == first(gtb.b[[4, 5]])
  @test ngtb.b[3] == first(gtb.b[[1, 2]])

  ngtb = gtb |> Aggregate(PointSet(pts1), :a => median, :b => last)
  @test domain(ngtb) == PointSet(pts1)
  @test ngtb.a[1] == gtb.a[3]
  @test ngtb.a[2] == median(gtb.a[[4, 5]])
  @test ngtb.a[3] == median(gtb.a[[1, 2]])
  @test ngtb.b[1] == gtb.b[3]
  @test ngtb.b[2] == last(gtb.b[[4, 5]])
  @test ngtb.b[3] == last(gtb.b[[1, 2]])
end
