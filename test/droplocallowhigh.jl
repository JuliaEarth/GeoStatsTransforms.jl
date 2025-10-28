@testset "DropLocalLowHigh" begin
  @test !isrevertible(DropLocalLowHigh(1))

  # basic tests
  gtb = georef((; a=[1.0, 2.0, 3.0, 100.0, 5.0]))
  ntb = gtb |> DropLocalLowHigh(1)
  @test ntb.a == [2.0, 3.0]
  ntb = gtb |> DropLocalLowHigh(1, low=0.0, high=0.98)
  @test ntb.a == [1.0, 2.0, 3.0, 5.0]
  ntb = gtb |> DropLocalLowHigh(1, low=0.02, high=1.0)
  @test ntb.a == [2.0, 3.0, 100.0]

  # tests with aliases
  gtb = georef((; a=[1.0, 2.0, 3.0, 100.0, 5.0]))
  ntb = gtb |> DropLocalMinima(1)
  @test ntb.a == [2.0, 3.0, 100.0]
  ntb = gtb |> DropLocalMaxima(1)
  @test ntb.a == [1.0, 2.0, 3.0, 5.0]

  # missing values
  gtb = georef((; a=[1.0, missing, 3.0, 100.0, 5.0]))
  ntb = gtb |> DropLocalLowHigh(1)
  @test isequal(ntb.a, [1.0, missing])
end
