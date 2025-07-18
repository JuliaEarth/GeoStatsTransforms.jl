@testset "Quenching" begin
  f = GaussianTransiogram()
  @test !isrevertible(Quenching(f))

  # basic test
  f = SphericalTransiogram()
  d = georef((; z=[iseven(i) for i in 1:10, j in 1:10]))
  q = d |> Quenching(f)
  @test Set(q.z) == Set([true, false])

  # skip specific locations
  skip = [1, 10, 10 * 10 - 10 + 1, 10 * 10]
  q = d |> Quenching(f; skip)
  @test all(q.z[skip] .== d.z[skip])
end
