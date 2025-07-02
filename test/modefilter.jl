@testset "ModeFilter" begin
  f = GaussianTransiogram()
  @test !isrevertible(ModeFilter())

  # basic test
  d = georef((; z=[iseven(i) for i in 1:10, j in 1:10]))
  p = d |> ModeFilter()
  @test Set(p.z) == Set([true, false])

  # skip specific locations
  skip = [1, 10, 10 * 10 - 10 + 1, 10 * 10]
  p = d |> ModeFilter(; skip)
  @test all(p.z[skip] .== d.z[skip])
end
