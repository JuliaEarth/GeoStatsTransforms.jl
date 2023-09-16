@testset "Feature" begin
  # transforms that are revertible
  d = georef((z=rand(100), w=rand(100)))
  for p in [
    Select(:z),
    Reject(:z),
    Rename(:z => :a),
    StdNames(),
    Sort(:z),
    Sample(10),
    Filter(x -> true),
    DropMissing(),
    Replace(1.0 => 2.0),
    Coalesce(value=0.0),
    Coerce(:z => Continuous),
    Identity(),
    Center(),
    Scale(),
    MinMax(),
    Interquartile(),
    ZScore(),
    Quantile(),
    Functional(sin),
    EigenAnalysis(:V),
    PCA(),
    DRS(),
    SDS(),
    RowTable(),
    ColTable()
  ]
    n, c = apply(p, d)
    t = Tables.columns(n)
    r = revert(p, n, c)
    @test n isa AbstractGeoTable
    @test r isa AbstractGeoTable
  end

  # transforms with categorical variables
  d = georef((c=categorical([1, 2, 3]),))
  for p in [Levels(:c => [1, 2, 3]), OneHot(:c)]
    n, c = apply(p, d)
    t = Tables.columns(n)
    r = revert(p, n, c)
    @test n isa AbstractGeoTable
    @test r isa AbstractGeoTable
  end

  d = georef((z=rand(100), w=rand(100)))
  p = Select(:w)
  n, c = apply(p, d)
  t = Tables.columns(n)
  @test Tables.columnnames(t) == (:w, :geometry)

  d = georef((z=rand(100), w=rand(100)))
  p = Sample(100)
  n, c = apply(p, d)
  r = revert(p, n, c)
  @test r == d
  t = Tables.columns(n)
  @test Tables.columnnames(t) == (:z, :w, :geometry)

  d = georef((a=[1, missing, 3], b=[3, 2, 1]))
  p = DropMissing()
  n, c = apply(p, d)
  @test Tables.columns(values(n)) == (a=[1, 3], b=[3, 1])
  @test nelements(domain(n)) == 2
  r = revert(p, n, c)
  @test r == d
end
