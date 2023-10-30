@testset "Simulate" begin
  a = rand(100)
  b = rand(100)
  c = rand(100)
  gtb = georef((; a, b, c), CartesianGrid(10, 10))

  dom = CartesianGrid(20, 20)
  sim = gtb |> Simulate(dom, 5, :a => GaussianProcess())
  @test nrow(sim) == 400
  @test propertynames(sim) == [:a_1, :a_2, :a_3, :a_4, :a_5, :geometry]
  sim = gtb |> Simulate(dom, 10, :a => GaussianProcess())
  @test nrow(sim) == 400
  @test propertynames(sim) == [:a_01, :a_02, :a_03, :a_04, :a_05, :a_06, :a_07, :a_08, :a_09, :a_10, :geometry]

  pts = rand(Point2, 200)
  sim = gtb |> Simulate(pts, 2, [:b, :c] => GaussianProcess())
  @test nrow(sim) == 200
  @test propertynames(sim) == [:b_1, :c_1, :b_2, :c_2, :geometry]
  sim = gtb |> Simulate(pts, [:b, :c] => GaussianProcess())
  @test nrow(sim) == 200
  @test propertynames(sim) == [:b_1, :c_1, :geometry]
end
