@testset "Clustering" begin
  @testset "SLIC" begin
    Z = [ones(10, 10) 2ones(10, 10); 3ones(10, 10) 4ones(10, 10)]
    𝒮 = georef((Z=Z,))
    C = 𝒮 |> SLIC(4, 1.0)
    @test C.label == vec(Z')

    𝒮 = georef((z=[√(i^2 + j^2) for i in 1:100, j in 1:100],))
    C = 𝒮 |> SLIC(50, 0.001)
    @test 50 ≤ length(unique(C.label)) ≤ 60

    # test SLIC with heterogeneous data
    Z = (a=rand(10), b=1:10, x=rand(10), y=rand(10))
    𝒮 = georef(Z, (:x, :y))
    C = 𝒮 |> SLIC(2, 1.0)
    @test domain(C) == domain(𝒮)
    @test Set(C.label) ⊆ Set(1:2)

    # test SLIC for orphaned points
    a = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    x = [
      0.4993029939801461,
      0.14954882636793432,
      0.23118957975519616,
      0.6816610871344635,
      0.6665309965318731,
      0.691522274292691,
      0.012495903053589608,
      0.9831177095525963,
      0.4445263730141056,
      0.2175871587746574
    ]
    y = [
      0.32721108209880256,
      0.11427387079564899,
      0.826401075107011,
      0.6164294766961782,
      0.6562529361193601,
      0.43388375115444644,
      0.7624847842129086,
      0.1516623758764959,
      0.07641616063237144,
      0.8669098569279463
    ]
    Z = (a=a, x=x, y=y)
    𝒮 = georef(Z, (:x, :y))
    C = 𝒮 |> SLIC(2, 1.0)
    @test Set(C.label) ⊆ Set(1:2)

    # test SLIC with weights in attribute columns
    z1 = [√((i - 0)^2 + (j - 0)^2) for i in 1:100, j in 1:100]
    z2 = [√((i - 100)^2 + (j - 100)^2) for i in 1:100, j in 1:100]
    𝒮 = georef((z1=z1, z2=z2))
    w1 = Dict(:z1 => 10, :z2 => 0.1)
    w2 = Dict(:z1 => 0.1, :z2 => 10)
    C1 = 𝒮 |> SLIC(50, 0.001, weights=w1)
    C2 = 𝒮 |> SLIC(50, 0.001, weights=w2)
    @test 50 ≤ length(unique(C1.label)) ≤ 60
    @test 50 ≤ length(unique(C2.label)) ≤ 60

    # test GeoClustering.slic_srecursion function
    k = 20
    l = [10.0, 100.0, 1000.0]
    s = GeoStatsTransforms.slic_srecursion(k, l)
    @test s[1] == 10 / 3 && s[2] == 100 / 3 && s[3] == 1000 / 3

    # the following test deals with the case where the bounding box
    # of the data has very different sides, one of which is too small
    # we want to make sure that the initialization of centroids always
    # returns a non-empty set
    k = 1
    m = 0.000001
    x = LinRange(550350.6224548942, 552307.2106300013, 1200)
    y = LinRange(9.35909841165263e6, 9.36050447440832e6, 1200)
    z = LinRange(-44.90690201082941, 351.4007207008662, 1200)
    𝒟 = PointSet(collect(zip(x, y, z)))
    s = GeoStatsTransforms.slic_spacing(𝒟, k)
    lo, up = to.(extrema(boundingbox(𝒟)))
    ranges = [(l + sᵢ / 2):sᵢ:u for (l, sᵢ, u) in zip(lo, s, up)]
    @test !isempty(Iterators.product(ranges...))
    c = GeoStatsTransforms.slic_initialization(𝒟, s)
    @test !isempty(c)
  end

  @testset "GHC" begin
    Z = [ones(10, 10) 2ones(10, 10); 3ones(10, 10) 4ones(10, 10)]
    𝒮 = georef((Z=Z,))
    C = 𝒮 |> GHC(4, 1.0)
    𝒮′ = georef(values(𝒮), centroid.(domain(𝒮)))
    C′ = 𝒮′ |> GHC(4, 1.0)
    @test C.label == categorical(vec(Z'))
    @test C.label == C′.label

    𝒮 = georef((z=[√(i^2 + j^2) for i in 1:10, j in 1:10],))
    C = 𝒮 |> GHC(10, 1.0)
    @test length(unique(C.label)) == 10

    𝒮 = georef((z=[10sin(i / 10) + j for i in 1:10, j in 1:10],))
    C = 𝒮 |> GHC(3, 1.0)
    @test length(unique(C.label)) == 3

    # multiple k
    𝒮 = georef((z=[10sin(i / 10) + j for i in 1:10, j in 1:10],))
    C = 𝒮 |> GHC([3, 5], 1.0)
    @test length(unique(C.label1)) == 3
    @test length(unique(C.label2)) == 5

    # objects in Hilbert space
    𝒮 = georef((z=rand(Composition{3}, 100),))
    C = 𝒮 |> GHC(5, 1.0)
    @test length(unique(C.label)) == 5
  end

  @testset "GSC" begin
    𝒮 = georef((Z=[10sin(i / 10) + j for i in 1:10, j in 1:10],))
    C = 𝒮 |> GSC(10, 2.0)
    @test Set(C.label) == Set(1:10)
  end
end
