@testset "MaxPosterior" begin
  f = GaussianTransiogram()
  @test !isrevertible(MaxPosterior(f))

  # basic test
  f = SphericalTransiogram()
  d = georef((; z=[iseven(i) for i in 1:10, j in 1:10]))
  p = d |> MaxPosterior(f)
  @test Set(p.z) == Set([true, false])

  # skip specific locations
  skip = [1, 10, 10 * 10 - 10 + 1, 10 * 10]
  p = d |> MaxPosterior(f; skip)
  @test all(p.z[skip] .== d.z[skip])
end
