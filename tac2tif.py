import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import calculate_default_transform, reproject, Resampling

tac_file = "RAF20.tac"
out_tif = "RAF20_epsg4326.tif"
out_tif_2154 = "RAF20_epsg2154.tif"

with open(tac_file) as f:
    header = f.readline().split()

    min_lon  = float(header[0])
    max_lon  = float(header[1])
    min_lat  = float(header[2])
    max_lat  = float(header[3])
    step_lon = float(header[4])
    step_lat = float(header[5])

    values = []

    for line in f:
        parts = line.split()
        for i in range(0, len(parts), 2):
            values.append(float(parts[i]))  # ignore flag

values = np.array(values)

ncols = int(round((max_lon - min_lon) / step_lon)) + 1
nrows = int(round((max_lat - min_lat) / step_lat)) + 1

grid = values.reshape((nrows, ncols))

transform = from_origin(
    min_lon,
    max_lat,
    step_lon,
    step_lat
)

with rasterio.open(
    out_tif,
    "w",
    driver="GTiff",
    height=nrows,
    width=ncols,
    count=1,
    dtype="float32",
    crs="EPSG:4326",
    transform=transform
) as dst:
    dst.write(grid.astype("float32"), 1)

print("RAF20 converti en GeoTIFF :", out_tif)

# Reprojection en EPSG:2154 avec interpolation bicubique
crs_2154 = "EPSG:2154"

with rasterio.open(out_tif) as src:
    transform_2154, width_2154, height_2154 = calculate_default_transform(
        src.crs, crs_2154, src.width, src.height, *src.bounds
    )
    profile = src.meta.copy()
    profile.update({
        "crs": crs_2154,
        "transform": transform_2154,
        "width": width_2154,
        "height": height_2154,
    })

    with rasterio.open(out_tif_2154, "w", **profile) as dst:
        reproject(
            source=rasterio.band(src, 1),
            destination=rasterio.band(dst, 1),
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform_2154,
            dst_crs=crs_2154,
            resampling=Resampling.cubic,
            src_nodata=0.0,
            dst_nodata=0.0
        )

print("RAF20 reprojeté en EPSG:2154 :", out_tif_2154)
