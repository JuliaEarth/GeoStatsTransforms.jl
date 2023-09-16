@testset "Rasterize" begin
  @test isrevertible(Rasterize(10, 10)) == true

  x = [1, 2, 3, 4, 5]
  y = [1.1, 2.2, 3.3, 4.4, 5.5]
  pts = Point2[(3, 9), (7, 8), (8, 5), (5, 4), (1, 5)]
  seg1 = Segment(pts[1], pts[2])
  seg2 = Segment(pts[2], pts[3])
  seg3 = Segment(pts[3], pts[4])
  seg4 = Segment(pts[4], pts[5])
  seg5 = Segment(pts[5], pts[1])
  poly1 = PolyArea((2, 0), (6, 2), (2, 2))
  poly2 = PolyArea((0, 6), (3, 8), (0, 10))
  poly3 = PolyArea((3, 6), (9, 6), (9, 9), (6, 9))
  poly4 = PolyArea((7, 0), (10, 0), (10, 4), (7, 4))
  poly5 = PolyArea((1, 3), (5, 3), (6, 6), (3, 8), (0, 6))

  gtb = georef((; x, y), pts)
  grid = CartesianGrid(10, 10)
  trans = Rasterize(grid)
  ngtb, cache = apply(trans, gtb)
  linds = LinearIndices((10, 10))
  @test ngtb.x[linds[3, 9]] == 1
  @test ngtb.x[linds[7, 8]] == 2
  @test ngtb.x[linds[8, 5]] == 3
  @test ngtb.x[linds[5, 4]] == 4
  @test ngtb.x[linds[1, 5]] == 5
  @test ngtb.y[linds[3, 9]] == 1.1
  @test ngtb.y[linds[7, 8]] == 2.2
  @test ngtb.y[linds[8, 5]] == 3.3
  @test ngtb.y[linds[5, 4]] == 4.4
  @test ngtb.y[linds[1, 5]] == 5.5

  gtb = georef((; x, y), [seg1, seg2, seg3, seg4, seg5])
  grid = CartesianGrid((0, 0), (10, 10), dims=(20, 20))
  trans = Rasterize(grid)
  ngtb, cache = apply(trans, gtb)
  linds = LinearIndices((20, 20))
  @test ngtb.x[linds[10, 17]] == 1
  @test ngtb.x[linds[15, 13]] == 2
  @test ngtb.x[linds[13, 9]] == 3
  @test ngtb.x[linds[6, 9]] == 4
  @test ngtb.x[linds[4, 14]] == 5
  @test ngtb.y[linds[10, 17]] == 1.1
  @test ngtb.y[linds[15, 13]] == 2.2
  @test ngtb.y[linds[13, 9]] == 3.3
  @test ngtb.y[linds[6, 9]] == 4.4
  @test ngtb.y[linds[4, 14]] == 5.5

  gtb = georef((; x, y), [poly1, poly2, poly3, poly4, poly5])
  trans = Rasterize(20, 20, :x => last, :y => std)
  ngtb, cache = apply(trans, gtb)
  linds = LinearIndices((20, 20))
  @test ngtb.x[linds[7, 3]] == 1
  @test ngtb.x[linds[3, 16]] == 2
  @test ngtb.x[linds[15, 15]] == 3
  @test ngtb.x[linds[17, 5]] == 4
  @test ngtb.x[linds[6, 11]] == 5
  @test ngtb.y[linds[7, 3]] == 1.1
  @test ngtb.y[linds[3, 16]] == 2.2
  @test ngtb.y[linds[15, 15]] == 3.3
  @test ngtb.y[linds[17, 5]] == 4.4
  @test ngtb.y[linds[6, 11]] == 5.5
  # intersection: poly3 with poly5
  @test ngtb.x[linds[9, 13]] == last(gtb.x[[3, 5]])
  @test ngtb.y[linds[9, 13]] == std(gtb.y[[3, 5]])

  # revert
  gtb = georef((; z=1:4), [poly1, poly2, poly3, poly4])
  trans = Rasterize(200, 200)
  ngtb, cache = apply(trans, gtb)
  rgtb = revert(trans, ngtb, cache)
  inds = filter(!iszero, unique(cache))
  @test isapprox(area(gtb.geometry[inds[1]]), area(rgtb.geometry[1]), atol=0.5)
  @test isapprox(area(gtb.geometry[inds[2]]), area(rgtb.geometry[2]), atol=0.5)
  @test isapprox(area(gtb.geometry[inds[3]]), area(rgtb.geometry[3]), atol=0.5)
  @test isapprox(area(gtb.geometry[inds[4]]), area(rgtb.geometry[4]), atol=0.5)
  # geotable with "mask" column
  gtb = georef((; z=1:4, mask=4:-1:1), [poly1, poly2, poly3, poly4])
  trans = Rasterize(10, 10)
  ngtb, cache = apply(trans, gtb)
  rgtb = revert(trans, ngtb, cache)
  @test nrow(rgtb) == nrow(gtb)
  @test ncol(rgtb) == ncol(gtb)
  @test propertynames(rgtb) == propertynames(gtb)
end
