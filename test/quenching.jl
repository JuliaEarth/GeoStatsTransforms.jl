@testset "Quenching" begin
  f = SphericalTransiogram()
  @test !isrevertible(Quenching(f))

  # basic test
  f = GaussianTransiogram()
  d = georef((; z=[iseven(i) for i in 1:(100 * 100)]), CartesianGrid(100, 100))
  q = d |> Quenching(f)
  @test Set(q.z) == Set([true, false])

  # skip specific locations
  skip = [1, 100, 100 * 100 - 100 + 1, 100 * 100]
  q = d |> Quenching(f; skip)
  @test all(q.z[skip] .== d.z[skip])
end
