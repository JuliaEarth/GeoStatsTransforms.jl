@testset "Potrace" begin
  # challenging case with letters
  img = load(joinpath(datadir, "letters.png"))
  gtb = georef((color=img,))
  trans = Potrace(1)
  ngtb, cache = apply(trans, gtb)
  ndom = domain(ngtb)
  @test nelements(ndom) == 2
  @test eltype(ndom) <: Multi
  polys1 = parent(ndom[1])
  polys2 = parent(ndom[2])
  @test length(polys1) == 4
  @test length(polys2) == 2
  rgtb = revert(trans, ngtb, cache)
  dom = domain(gtb)
  rdom = domain(rgtb)
  @test rdom isa Grid
  @test size(rdom) == size(dom)
  @test minimum(rdom) == minimum(dom)
  @test maximum(rdom) == maximum(dom)
  @test spacing(rdom) == spacing(dom)

  # concentric circles
  ball1 = Ball((0, 0), 1)
  ball2 = Ball((0, 0), 2)
  ball3 = Ball((0, 0), 3)
  grid = CartesianGrid((-5, -5), (5, 5), dims=(100, 100))
  inds1 = centroid.(grid) .∈ Ref(ball1)
  inds2 = centroid.(grid) .∈ Ref(ball2)
  inds3 = centroid.(grid) .∈ Ref(ball3)
  mask = zeros(100, 100)
  mask[inds3] .= 1
  mask[inds2] .= 0
  mask[inds1] .= 1
  dat = georef((mask=mask,))
  new = dat |> Potrace(1)
  dom = domain(new)
  @test nelements(dom) == 2
  @test eltype(dom) <: Multi
  polys1 = parent(dom[1])
  polys2 = parent(dom[2])
  @test length(polys1) == 2
  @test length(polys2) == 2
  new1 = dat |> Potrace(1, ϵ=0.1)
  new2 = dat |> Potrace(1, ϵ=0.5)
  dom1 = domain(new1)
  dom2 = domain(new2)
  for (g1, g2) in zip(dom1, dom2)
    @test nvertices(g1) > nvertices(g2)
  end

  # make sure that aggregation works
  Z = [sin(i / 10) + sin(j / 10) for i in 1:100, j in 1:100]
  M = Z .> 0
  Ω = georef((Z=Z, M=M))
  𝒯 = Ω |> Potrace(:M, :Z => mean)
  masks = unique(Ω.M)
  @test nelements(domain(𝒯)) == 2
  @test Set(𝒯.M) == Set([true, false])
  @test all(z -> -1 ≤ z ≤ 1, 𝒯.Z)
  @test 𝒯.Z[1] == mean(Ω.Z[masks[1] .== Ω.M])
  @test 𝒯.Z[2] == mean(Ω.Z[masks[2] .== Ω.M])

  # units
  gtb = georef((; T=Z * u"K", M))
  ngtb = gtb |> Potrace(:M)
  masks = unique(gtb.M)
  @test GeoStatsTransforms.elunit(ngtb.T) == u"K"
  @test ngtb.T[1] ≈ mean(gtb.T[masks[1] .== gtb.M])
  @test ngtb.T[2] ≈ mean(gtb.T[masks[2] .== gtb.M])

  # affine units
  gtb = georef((; T=Z * u"°C", M))
  ngtb = gtb |> Potrace(:M)
  masks = unique(gtb.M)
  @test GeoStatsTransforms.elunit(ngtb.T) == u"K"
  v = GeoStatsTransforms.uadjust(gtb.T[masks[1] .== gtb.M])
  @test ngtb.T[1] ≈ mean(v)
  v = GeoStatsTransforms.uadjust(gtb.T[masks[2] .== gtb.M])
  @test ngtb.T[2] ≈ mean(v)
end
