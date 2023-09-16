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
  x = rand(100)
  y = rand(1:10, 100)
  table = (; x, y)
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
    @test ndata.x[i] == mean(sdata.x[(j - 9):j])
  end

  for i in 1:10
    j = i * 10
    @test ndata.y[i] == first(sdata.y[(j - 9):j])
  end

  # custom aggregators
  # colspec: index
  ndata = sdata |> UniqueCoords(1 => std, 2 => median)
  @test nrow(ndata) == 10

  for i in 1:10
    j = i * 10
    @test ndata.x[i] == std(sdata.x[(j - 9):j])
  end

  for i in 1:10
    j = i * 10
    @test ndata.y[i] == median(sdata.y[(j - 9):j])
  end

  # colspec: symbols
  ndata = sdata |> UniqueCoords(:x => last, :y => first)
  @test nrow(ndata) == 10

  for i in 1:10
    j = i * 10
    @test ndata.x[i] == last(sdata.x[(j - 9):j])
  end

  for i in 1:10
    j = i * 10
    @test ndata.y[i] == first(sdata.y[(j - 9):j])
  end

  # colspec: strings
  ndata = sdata |> UniqueCoords("x" => maximum, "y" => minimum)
  @test nrow(ndata) == 10

  for i in 1:10
    j = i * 10
    @test ndata.x[i] == maximum(sdata.x[(j - 9):j])
  end

  for i in 1:10
    j = i * 10
    @test ndata.y[i] == minimum(sdata.y[(j - 9):j])
  end
end
