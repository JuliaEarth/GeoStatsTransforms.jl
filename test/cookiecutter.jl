@testset "CookieCutter" begin
  table = (facies=[1, 0, 1], poro=[0.5, 0.9, 0.1])
  coord = [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)]
  geotable = georef(table, coord)

  trainimg = geostatsimage("Strebelle")
  parent = QuiltingProcess(trainimg=trainimg, tilesize=(30, 30))

  child0 = GaussianProcess(variogram=SphericalVariogram(range=20.0, sill=0.2))
  child1 = GaussianProcess(variogram=SphericalVariogram(MetricBall((200.0, 20.0))))

  sdomain = CartesianGrid(100, 100)
  transform = CookieCutter(sdomain, :facies => parent, :poro => [0 => child0, 1 => child1])
  @test !isrevertible(transform)
  ngtb = transform(geotable)
  @test nrow(ngtb) == 10000
  @test propertynames(ngtb) == [:facies_1, :poro_1, :geometry]

  transform = CookieCutter(sdomain, 3, :facies => parent, :poro => [0 => child0, 1 => child1])
  @test !isrevertible(transform)
  ngtb = transform(geotable)
  @test nrow(ngtb) == 10000
  @test propertynames(ngtb) == [:facies_1, :facies_2, :facies_3, :poro_1, :poro_2, :poro_3, :geometry]

  # throw: CookieCutter without children
  @test_throws ArgumentError CookieCutter(sdomain, :facies => parent)
  @test_throws ArgumentError CookieCutter(sdomain, 3, :facies => parent)
  # throw: invalid child map
  @test_throws ArgumentError CookieCutter(sdomain, :facies => parent, :poro => [0 => nothing, 1 => child1])
  @test_throws ArgumentError CookieCutter(sdomain, 3, :facies => parent, :poro => [0 => nothing, 1 => child1])
end
