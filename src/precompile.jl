# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using PrecompileTools

@setup_workload begin
  gtb = georef((; Z=rand(10, 10)))
  proc = GaussianProcess()
  @compile_workload begin
    gtb |> Interpolate(gtb.geometry)
    gtb |> InterpolateNeighbors(gtb.geometry)
    gtb |> InterpolateMissing()
    gtb |> InterpolateNaN()
    gtb |> Simulate(gtb.geometry, :Z => proc)
    gtb |> UniqueCoords()
    gtb |> Aggregate(gtb.geometry)
    gtb |> Transfer(gtb.geometry)
    gtb |> Upscale(2, 2)
    gtb |> Downscale(2, 2)
    gtb |> SLIC(3, 1.0, as=:cluster)
    gtb |> GHC(3, 1.0, as=:cluster)
    gtb |> GSC(3, 2.0, as=:cluster)
    gtb |> Rasterize(gtb.geometry)
    gtb |> Potrace(:Z)
    gtb |> Detrend(:Z)
  end
end
