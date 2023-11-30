@testset "UniqueCoords" begin
  @test isrevertible(UniqueCoords()) == false

  X = [i * j for i in 1:2, j in 1:1_000_000]
  z = rand(1_000_000)
  d = georef((z=[z; z],), [X X])
  u = d |> UniqueCoords()
  du = domain(u)
  p = [centroid(du, i) for i in 1:nelements(du)]
  U = reduce(hcat, coordinates.(p))
  @test nelements(du) == 1_000_000
  @test Set(eachcol(U)) == Set(eachcol(X))

  X = rand(3, 100)
  z = rand(100)
  n = [string(i) for i in 1:100]
  Xd = hcat(X, X[:, 1:10])
  zd = vcat(z, z[1:10])
  nd = vcat(n, n[1:10])
  sdata = georef((z=zd, n=nd), PointSet(Xd))
  ndata = sdata |> UniqueCoords()
  @test nrow(ndata) == 100

  # domain with repeated points
  a = rand(100)
  b = rand(1:10, 100)
  table = (; a, b)
  points = Point.(1:10, 10:-1:1)

  pset = PointSet(rand(points, 100))
  sdata = georef(table, pset)
  ndata = sdata |> UniqueCoords()
  @test nrow(ndata) == 10

  # aggregators
  pset = PointSet(repeat(points, inner=10))
  sdata = georef(table, pset)

  # default aggregators
  ndata = sdata |> UniqueCoords()
  @test nrow(ndata) == 10

  for i in 1:10
    j = i * 10
    @test ndata.a[i] == mean(sdata.a[(j - 9):j])
  end

  for i in 1:10
    j = i * 10
    @test ndata.b[i] == first(sdata.b[(j - 9):j])
  end

  # custom aggregators
  # selector: indices
  ndata = sdata |> UniqueCoords(1 => std, 2 => median)
  @test nrow(ndata) == 10

  for i in 1:10
    j = i * 10
    @test ndata.a[i] == std(sdata.a[(j - 9):j])
  end

  for i in 1:10
    j = i * 10
    @test ndata.b[i] == median(sdata.b[(j - 9):j])
  end

  # selector: symbols
  ndata = sdata |> UniqueCoords(:a => last, :b => first)
  @test nrow(ndata) == 10

  for i in 1:10
    j = i * 10
    @test ndata.a[i] == last(sdata.a[(j - 9):j])
  end

  for i in 1:10
    j = i * 10
    @test ndata.b[i] == first(sdata.b[(j - 9):j])
  end

  # selector: strings
  ndata = sdata |> UniqueCoords("a" => maximum, "b" => minimum)
  @test nrow(ndata) == 10

  for i in 1:10
    j = i * 10
    @test ndata.a[i] == maximum(sdata.a[(j - 9):j])
  end

  for i in 1:10
    j = i * 10
    @test ndata.b[i] == minimum(sdata.b[(j - 9):j])
  end

  # units
  sdata = georef((; T=rand(100) * u"K"), pset)
  ndata = sdata |> UniqueCoords()
  @test nrow(ndata) == 10
  @test unit(eltype(ndata.T)) == u"K"

  for i in 1:10
    j = i * 10
    @test ndata.T[i] == mean(sdata.T[(j - 9):j])
  end

  # affine units
  sdata = georef((; T=rand(100) * u"°C"), pset)
  ndata = sdata |> UniqueCoords()
  @test nrow(ndata) == 10
  @test unit(eltype(ndata.T)) == u"K"

  for i in 1:10
    j = i * 10
    v = GeoStatsTransforms._absunit(sdata.T[(j - 9):j])
    @test ndata.T[i] == mean(v)
  end

  # units and missings
  sdata = georef((; T=shuffle([fill(missing, 50); rand(50)]) * u"K"), pset)
  ndata = sdata |> UniqueCoords()
  @test nrow(ndata) == 10
  @test unit(eltype(ndata.T)) == u"K"

  for i in 1:10
    j = i * 10
    v = GeoStatsTransforms._skipmissing(mean)(sdata.T[(j - 9):j])
    @test isequal(ndata.T[i], v)
  end
end
