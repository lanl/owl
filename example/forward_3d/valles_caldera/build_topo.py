
import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import elevation
import pyproj
import matplotlib.pyplot as plt
import os
import math

from rasterio.warp import calculate_default_transform
from rasterio.transform import xy

def write_array(w, filename, dtype=np.float32, ascii=False):

    w = np.asarray(w, dtype=dtype)
    w = np.transpose(w)

    if ascii:
        np.savetxt(filename, w.transpose(), fmt='%e', delimiter=' ')
        # w.tofile(filename, sep='\n', format='%e')
    else:
        w.tofile(filename)

# -------------------------
# 1. User parameters
# -------------------------
# lat_start, lat_end = 35.788235, 36.056627   # degrees

lat_start, lat_end = 35.75, 36.05   # degrees
lon_start, lon_end = -106.738299, -106.270029  # degrees

resolution_m = 50                 # desired resolution in meters (projected)

output_dem = 'dem.tif'  # final DEM in current directory

# -------------------------
# 2. Configure elevation to use current directory
# -------------------------
cache_dir = os.path.join(os.getcwd(), 'srtm_cache')
os.makedirs(cache_dir, exist_ok=True)

# -------------------------
# 3. Download and clip SRTM to current directory
# -------------------------
elevation.clip(bounds=(lon_start, lat_start, lon_end, lat_end),
               output=output_dem, product='SRTM1', cache_dir=cache_dir)

# -------------------------
# 4. Open DEM
# -------------------------
with rasterio.open(cache_dir + '/SRTM1/' + output_dem) as src:
    src_crs = src.crs
    src_transform = src.transform
    src_width = src.width
    src_height = src.height

    print(src_width, src_height)

    # -------------------------
    # 5. Estimate lat/lon resolution corresponding to target meters
    # -------------------------
    center_lat = (lat_start + lat_end) / 2
    meters_per_deg_lat = 111132
    meters_per_deg_lon = 111320 * math.cos(math.radians(center_lat))

    res_deg_lat = resolution_m / meters_per_deg_lat
    res_deg_lon = resolution_m / meters_per_deg_lon

    new_width = int((lon_end - lon_start) / res_deg_lon)
    new_height = int((lat_end - lat_start) / res_deg_lat)

    print(new_width, new_height)

    # -------------------------
    # 6. Resample DEM in lat/lon
    # -------------------------
    dem_resampled_latlon = src.read(
        1,
        out_shape=(new_height, new_width),
        resampling=rasterio.enums.Resampling.bilinear
    )

    new_transform = rasterio.transform.from_bounds(lon_start, lat_start, lon_end, lat_end,
                                                   new_width, new_height)

# -------------------------
# 7. Reproject to UTM (meters)
# -------------------------
center_lon = (lon_start + lon_end) / 2
utm_zone = int(center_lon // 6 + 31)
print(utm_zone)
proj_utm = pyproj.CRS(f"+proj=utm +zone={utm_zone} +datum=WGS84 +units=m +no_defs")

# Prepare destination array
dem_utm = np.empty_like(dem_resampled_latlon)

reproject(
    source=dem_resampled_latlon,
    destination=dem_utm,
    src_transform=new_transform,
    src_crs=src_crs,
    dst_transform=None,   # rasterio will compute automatically
    dst_crs=proj_utm,
    resampling=Resampling.bilinear
)

# -------------------------
# 8. Plot result
# -------------------------
plt.figure(figsize=(10, 8))
plt.imshow(dem_utm, origin='lower', cmap='terrain', vmin=0, vmax=3000)
plt.colorbar(label='Elevation (m)')
plt.title('Topography Map (UTM meters)')
plt.show()


print(dem_utm.shape)
write_array(dem_utm, './model/dem_utm.bin')

dst_transform, width, height = calculate_default_transform(
    src_crs, proj_utm, dem_resampled_latlon.shape[1], dem_resampled_latlon.shape[0], 
    *[lon_start, lat_start, lon_end, lat_end]
)

# Top-left corner
tl_x, tl_y = xy(dst_transform, 0, 0)
# Top-right
tr_x, tr_y = xy(dst_transform, 0, width-1)
# Bottom-left
bl_x, bl_y = xy(dst_transform, height-1, 0)
# Bottom-right
br_x, br_y = xy(dst_transform, height-1, width-1)

print("UTM corners (meters):")
print("Top-left:", tl_x, tl_y)
print("Top-right:", tr_x, tr_y)
print("Bottom-left:", bl_x, bl_y)
print("Bottom-right:", br_x, br_y)

# # The result is: 
# (666, 844)
# UTM corners (meters):
# Top-left: 342855.56048806757 3990867.2228574837
# Top-right: 385588.09589952044 3990867.2228574837
# Bottom-left: 342855.56048806757 3956991.4307338847
# Bottom-right: 385588.09589952044 3956991.4307338847
