#
# The script is generate a planar fault ensemble source for 
# testing OWL 3D simulation with complex source. 
#

import numpy as np
import pandas as pd
import pyvista as pv

def generate_fault_source(
    model_size_km=(10.0, 10.0, 4.0),
    fault_center_km=(5.0, 5.0, 2.0),
    fault_length_km=6.0,
    fault_width_km=3.0,
    strike_deg=45.0,
    dip_deg=60.0,
    dx_km=0.1,
    rupture_velocity_km_s=2.5,
    hypocenter_local_km=(0.0, 0.0),
    t0_shift_s=0.0,
    add_random_time_perturbation=False,
    random_time_std_s=0.02,
    seed=1234,
):
    """
    Generate a kinematic rupture source as point ensemble: x, y, z, t0.

    Coordinates:
        x, y, z are in km
        z is positive downward
        t0 is in seconds

    Fault geometry:
        strike direction is horizontal
        dip direction is perpendicular to strike and downward
        local fault coordinates are:
            s = along-strike coordinate
            d = down-dip coordinate

    hypocenter_local_km:
        (s_hypo, d_hypo) in km relative to the fault center.
    """

    rng = np.random.default_rng(seed)

    Lx, Ly, Lz = model_size_km
    cx, cy, cz = fault_center_km

    strike = np.deg2rad(strike_deg)
    dip = np.deg2rad(dip_deg)

    # Unit vector along strike
    strike_vec = np.array([
        np.sin(strike),
        np.cos(strike),
        0.0
    ])

    # Horizontal unit vector perpendicular to strike
    dip_horizontal_vec = np.array([
        np.cos(strike),
        -np.sin(strike),
        0.0
    ])

    # Down-dip unit vector, with z positive downward
    dip_vec = np.cos(dip) * dip_horizontal_vec + np.sin(dip) * np.array([0.0, 0.0, 1.0])

    # Local coordinates on the fault plane
    s_vals = np.arange(-fault_length_km / 2, fault_length_km / 2 + dx_km, dx_km)
    d_vals = np.arange(-fault_width_km / 2, fault_width_km / 2 + dx_km, dx_km)

    S, D = np.meshgrid(s_vals, d_vals, indexing="xy")

    points = (
        np.array([cx, cy, cz])[None, None, :]
        + S[:, :, None] * strike_vec[None, None, :]
        + D[:, :, None] * dip_vec[None, None, :]
    )

    x = points[:, :, 0].ravel()
    y = points[:, :, 1].ravel()
    z = points[:, :, 2].ravel()
    s = S.ravel()
    d = D.ravel()

    # Keep only points inside the model
    inside = (
        (x >= 0.0) & (x <= Lx) &
        (y >= 0.0) & (y <= Ly) &
        (z >= 0.0) & (z <= Lz)
    )

    x = x[inside]
    y = y[inside]
    z = z[inside]
    s = s[inside]
    d = d[inside]

    # Rupture starts from local hypocenter
    s0, d0 = hypocenter_local_km

    rupture_distance_km = np.sqrt((s - s0) ** 2 + (d - d0) ** 2)

    t0 = rupture_distance_km / rupture_velocity_km_s + t0_shift_s

    if add_random_time_perturbation:
        t0 += rng.normal(0.0, random_time_std_s, size=t0.shape)
        t0 = np.maximum(t0, 0.0)

    source = pd.DataFrame({
        "x_km": x,
        "y_km": y,
        "z_km": z,
        "t0_s": t0,
        "s_local_km": s,
        "d_local_km": d,
    })

    source = source.sort_values("t0_s").reset_index(drop=True)

    return source

def double_couple_moment_tensor_enu_down(
    strike_deg,
    dip_deg,
    rake_deg,
    M0=1.0,
):
    """
    Return double-couple moment tensor for coordinates:

        x = east
        y = north
        z = down

    Strike is clockwise from north.
    Dip is downward from horizontal.
    Rake is slip direction in the fault plane.

    Returns:
        3x3 moment tensor in x-y-z coordinates.
    """

    strike = np.deg2rad(strike_deg)
    dip = np.deg2rad(dip_deg)
    rake = np.deg2rad(rake_deg)

    # Unit vector along strike
    e_strike = np.array([
        np.sin(strike),
        np.cos(strike),
        0.0
    ])

    # Unit vector down dip
    e_dip = np.array([
        np.cos(strike) * np.cos(dip),
        -np.sin(strike) * np.cos(dip),
        np.sin(dip)
    ])

    # Fault normal; chosen so that e_strike x e_dip = normal
    n = np.cross(e_strike, e_dip)
    n = n / np.linalg.norm(n)

    # Slip direction from rake
    slip = np.cos(rake) * e_strike + np.sin(rake) * e_dip
    slip = slip / np.linalg.norm(slip)

    # Double-couple moment tensor
    M = M0 * (
        np.outer(slip, n) + np.outer(n, slip)
    )

    return M, slip, n

def visualize_fault_source_depth_up_pyvista(
    source,
    point_size=8,
    show_model_box=True,
    model_size_km=(10.0, 10.0, 4.0),
):
    """
    Visualize source with depth plotted downward visually.

    Original source convention:
        z_km > 0 means deeper

    Visualization convention:
        plot_z = -z_km, so deeper points appear lower.
    """

    points = source[["x_km", "y_km", "z_km"]].to_numpy().copy()
    points[:, 2] *= -1.0

    t0 = source["t0_s"].to_numpy()

    cloud = pv.PolyData(points)
    cloud["t0_s"] = t0

    p = pv.Plotter()

    p.add_mesh(
        cloud,
        scalars="t0_s",
        point_size=point_size,
        render_points_as_spheres=True,
        cmap="viridis",
        scalar_bar_args={"title": "t0 [s]"},
    )

    if show_model_box:
        Lx, Ly, Lz = model_size_km

        box = pv.Box(bounds=(0, Lx, 0, Ly, -Lz, 0))

        p.add_mesh(
            box,
            style="wireframe",
            color="black",
            line_width=1,
            opacity=0.25,
        )

    p.add_axes()
    p.show_grid(
        xtitle="x [km]",
        ytitle="y [km]",
        ztitle="-depth [km]",
    )

    # p.view_isometric()
    p.show()
    p.screenshot("fault.png") 

###################################################################################################

# Example usage
source = generate_fault_source(
    model_size_km=(10.0, 10.0, 4.0),
    fault_center_km=(5.0, 5.0, 2.0),
    fault_length_km=6.0,
    fault_width_km=3.0,
    strike_deg=45.0,
    dip_deg=60.0,
    dx_km=0.1,
    rupture_velocity_km_s=2.5,
    hypocenter_local_km=(-2.0, 0.0),
)

# source.to_csv("fault_source_x_y_z_t0.csv", index=False)
xyzt0 = source[["x_km", "y_km", "z_km", "t0_s"]].to_numpy()

np.savetxt(
    "fault_source_xyzt0.txt",
    xyzt0,
    fmt="%.6f %.6f %.6f %.6f"
)

print(source.head())
print()
print(f"Number of source points: {len(source)}")
print(f"Minimum t0: {source['t0_s'].min():.3f} s")
print(f"Maximum t0: {source['t0_s'].max():.3f} s")

visualize_fault_source_depth_up_pyvista(source)

# MT
M, slip, normal = double_couple_moment_tensor_enu_down(
    strike_deg=45.0,
    dip_deg=60.0,
    rake_deg=90.0,
    M0=1.0,
)

print("Moment tensor:")
print(M)

print("Slip vector:", slip)
print("Normal vector:", normal)