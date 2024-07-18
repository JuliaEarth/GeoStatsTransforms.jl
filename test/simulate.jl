@testset "Simulate" begin
  @test !isrevertible(Simulate(CartesianGrid(10, 10), :a => GaussianProcess()))

  a = rand(1000)
  b = rand(1000)
  c = rand(1000)
  gtb = georef((; a, b, c), CartesianGrid(10, 10, 10))

  dom = CartesianGrid(20, 20, 20)
  sim = gtb |> Simulate(dom, 5, :a => GaussianProcess())
  @test nrow(sim) == 8000
  @test propertynames(sim) == [:a_1, :a_2, :a_3, :a_4, :a_5, :geometry]
  sim = gtb |> Simulate(dom, 10, :a => GaussianProcess())
  @test nrow(sim) == 8000
  @test propertynames(sim) == [:a_01, :a_02, :a_03, :a_04, :a_05, :a_06, :a_07, :a_08, :a_09, :a_10, :geometry]

  pts = rand(Point, 200)
  sim = gtb |> Simulate(pts, 3, [:b, :c] => GaussianProcess())
  @test nrow(sim) == 200
  @test propertynames(sim) == [:b_1, :c_1, :b_2, :c_2, :b_3, :c_3, :geometry]
  sim = gtb |> Simulate(pts, [:b, :c] => GaussianProcess())
  @test nrow(sim) == 200
  @test propertynames(sim) == [:b_1, :c_1, :geometry]
end
