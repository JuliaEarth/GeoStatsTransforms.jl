@testset "CookieCutter" begin
  table = (facies=[1,0,1], poro=[0.5, 0.9, 0.1])
  coord = [(25.,25.), (50.,75.), (75.,50.)]
  geotable = georef(table, coord)

  trainimg = geostatsimage("Strebelle")
  mprocess = QuiltingProcess(trainimg=trainimg, tilesize=(30, 30))

  cprocess0 = GaussianProcess(variogram=SphericalVariogram(range=20., sill=.2))
  cprocess1 = GaussianProcess(variogram=SphericalVariogram(MetricBall((200.,20.))))

  sdomain = CartesianGrid(100, 100)
  transform = CookieCutter(sdomain, :facies => mprocess, :poro => [0 => cprocess0, 1 => cprocess1])
  @test !isrevertible(transform)
  ngtb = transform(geotable)
  @test nrow(ngtb) == 10000
  @test propertynames(ngtb) == [:facies_1, :poro_1, :geometry]

  transform = CookieCutter(sdomain, 3, :facies => mprocess, :poro => [0 => cprocess0, 1 => cprocess1])
  @test !isrevertible(transform)
  ngtb = transform(geotable)
  @test nrow(ngtb) == 10000
  @test propertynames(ngtb) == [:facies_1, :facies_2, :facies_3, :poro_1, :poro_2, :poro_3, :geometry]

  # throw: CookieCutter without children
  @test_throws ArgumentError CookieCutter(sdomain, :facies => mprocess)
  @test_throws ArgumentError CookieCutter(sdomain, 3, :facies => mprocess)
  # throw: invalid child map
  @test_throws ArgumentError CookieCutter(sdomain, :facies => mprocess, :poro => [0 => nothing, 1 => cprocess1])
  @test_throws ArgumentError CookieCutter(sdomain, 3, :facies => mprocess, :poro => [0 => nothing, 1 => cprocess1])
end
