@testset "Gradient" begin
  @test !isrevertible(Gradient(1))

  # constant field has zero gradient everywhere
  gtb = georef((; a=ones(5, 5)))
  grad = gtb |> Gradient(1)
  @test names(grad) == ["a_x", "a_y", "geometry"]
  @test unit(eltype(grad.a_x)) == u"m^-1"
  @test unit(eltype(grad.a_y)) == u"m^-1"
  @test all(==(0u"m^-1"), grad.a_x)
  @test all(==(0u"m^-1"), grad.a_y)

  # units are treated correctly
  gtb = georef((; a=ones(5, 5)u"K"))
  grad = gtb |> Gradient("a")
  @test unit(eltype(grad.a_x)) == u"K * m^-1"
  @test unit(eltype(grad.a_y)) == u"K * m^-1"
  @test all(==(0u"K * m^-1"), grad.a_x)
  @test all(==(0u"K * m^-1"), grad.a_y)

  # views of regular grids
  gtb = georef((; a=ones(5, 5)))
  vtb = view(gtb, 1:5)
  grad = vtb |> Gradient(1)
  @test all(==(0u"m^-1"), grad.a_x)
  @test all(==(-0.5u"m^-1"), grad.a_y)
  vtb = view(gtb, 21:25)
  grad = vtb |> Gradient(1)
  @test all(==(0u"m^-1"), grad.a_x)
  @test all(==(0.5u"m^-1"), grad.a_y)
end
