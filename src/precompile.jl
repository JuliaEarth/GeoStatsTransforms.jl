# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using PrecompileTools

@setup_workload begin
  gtb = georef((; z=rand(10, 10)))
  @compile_workload begin
    gtb |> Detrend(:z)
    gtb |> Upscale(2, 2)
    gtb |> Downscale(2, 2)
    gtb |> Interpolate(gtb.geometry)
    gtb |> InterpolateNeighbors(gtb.geometry)
    gtb |> Rasterize(gtb.geometry)
    gtb |> Potrace(:z)
  end
end
