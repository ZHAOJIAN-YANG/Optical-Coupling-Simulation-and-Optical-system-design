import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import os

wl_design = 10.0
ww_design = 50.0

dw_l = 1.0
dw_w = 1.0
dw_r = 1.0

dr_rl = 0.5
dr_rw = 0.5

dz_rl_plus = 2.0
dz_rl_minus = 2.0

alpha_deg = 30.0
alpha = np.deg2rad(alpha_deg)

reflector_depth_ratio = 0.90

profile_points = 800
phi_points = 240

grid_points = 201
EXPORT_GRD = True

wr_design_0 = (
    ww_design
    - dw_r
    - dw_w
    - 2.0 * dr_rw
)

if wr_design_0 <= 0:
    raise ValueError(
        "wr_design_0 <= 0，請檢查波導尺寸、公差與對位誤差。"
    )

rr_design_0 = wr_design_0 / 2.0

numerator_hrl = (
    wl_design
    + ww_design
    + dw_l
    - dw_w
    - 2.0 * dr_rw
    + 2.0 * dr_rl
)

denominator_hrl = 2.0 * np.tan(alpha)

hrl_design = (
    numerator_hrl / denominator_hrl
    + dz_rl_plus
)

if hrl_design <= 0:
    raise ValueError(
        "計算得到的 h_rl_design <= 0，參數不合理。"
    )

h_r = reflector_depth_ratio * hrl_design

if h_r >= hrl_design:
    raise ValueError(
        "h_r 必須小於 h_rl_design。"
    )


def Sr_lower(z, r):

    wr_z = 2.0 * r

    numerator = (
        wl_design
        + wr_z
        + dw_l
        + dw_r
        + 2.0 * dr_rl
    )

    denominator = (
        2.0
        * (
            hrl_design
            + z
            - dz_rl_plus
        )
    )

    angle = np.arctan2(
        numerator,
        denominator
    )

    return np.tan(
        (angle - alpha) / 2.0
    )


def Sr_upper(z, r):

    wr_z = 2.0 * r

    numerator = (
        wl_design
        + wr_z
        - dw_l
        - dw_r
        - 2.0 * dr_rl
    )

    denominator = (
        2.0
        * (
            hrl_design
            + z
            - dz_rl_minus
        )
    )

    angle = np.arctan2(
        numerator,
        denominator
    )

    return np.tan(
        (angle + alpha) / 2.0
    )


def Sr_design(z, r):

    r_scalar = float(
        np.atleast_1d(r)[0]
    )

    wr_z = 2.0 * r_scalar

    numerator = (
        wl_design
        + wr_z
    )

    denominator = (
        2.0
        * (
            hrl_design
            + z
        )
    )

    angle = np.arctan2(
        numerator,
        denominator
    )

    return np.tan(
        0.5 * angle
    )


def check_step6_condition(z, r):

    r_scalar = float(
        np.atleast_1d(r)[0]
    )

    wr_z = 2.0 * r_scalar

    theta_nominal = np.arctan2(
        wl_design + wr_z,
        2.0 * (hrl_design + z)
    )

    theta_plus = np.arctan2(
        (
            wl_design
            + wr_z
            + dw_l
            + dw_r
            + 2.0 * dr_rl
        ),
        (
            2.0
            * (
                hrl_design
                + z
                - dz_rl_plus
            )
        )
    )

    theta_minus = np.arctan2(
        (
            wl_design
            + wr_z
            - dw_l
            - dw_r
            - 2.0 * dr_rl
        ),
        (
            2.0
            * (
                hrl_design
                + z
                + dz_rl_minus
            )
        )
    )

    delta_plus = (
        theta_plus
        - theta_nominal
    )

    delta_minus = (
        theta_nominal
        - theta_minus
    )

    condition_value = max(
        delta_plus,
        delta_minus
    )

    return (
        condition_value <= alpha,
        condition_value
    )


def reflector_ode(z, r):

    return [
        Sr_design(z, r)
    ]


z_eval = np.linspace(
    0.0,
    -h_r,
    profile_points
)

solution = solve_ivp(
    reflector_ode,
    (0.0, -h_r),
    [rr_design_0],
    t_eval=z_eval,
    method="RK45",
    rtol=1e-9,
    atol=1e-11
)

if not solution.success:
    raise RuntimeError(
        solution.message
    )

z_profile = solution.t
r_profile = solution.y[0]

z_plot = z_profile[::-1]
r_plot = r_profile[::-1]

Z_display = (
    z_plot
    + h_r
)

Sr_values = np.array([
    Sr_design(z, r)
    for z, r in zip(
        z_profile,
        r_profile
    )
])

Sr_min_values = np.array([
    Sr_lower(z, r)
    for z, r in zip(
        z_profile,
        r_profile
    )
])

Sr_max_values = np.array([
    Sr_upper(z, r)
    for z, r in zip(
        z_profile,
        r_profile
    )
])

condition_results = [
    check_step6_condition(z, r)
    for z, r in zip(
        z_profile,
        r_profile
    )
]

condition_ok = np.array(
    [x[0] for x in condition_results]
)

condition_value = np.array(
    [x[1] for x in condition_results]
)

wr_profile = (
    2.0
    * r_profile
)

integral_Sr = (
    r_profile[-1]
    - rr_design_0
)

step9_rhs = (
    2.0 * integral_Sr
    + wr_design_0
    - dw_l
    - dw_r
    - 2.0 * dr_rl
)

step9_ok = (
    wl_design
    <=
    step9_rhs
)

phi = np.linspace(
    0.0,
    2.0 * np.pi,
    phi_points
)

PHI, Z3D = np.meshgrid(
    phi,
    Z_display
)

R3D = np.tile(
    r_plot[:, None],
    (1, phi_points)
)

X3D = (
    R3D
    * np.sin(PHI)
)

Y3D = (
    R3D
    * np.cos(PHI)
)

print("=" * 65)
print("碗型反射鏡面設計")
print("=" * 65)
print(
    f"Step 2  wr_design(0) = "
    f"{wr_design_0:.4f} μm"
)
print(
    f"         rr_design(0) = "
    f"{rr_design_0:.4f} μm"
)
print(
    f"Step 4  h_rl_design   = "
    f"{hrl_design:.4f} μm"
)
print(
    f"繪圖 reflector 深度 h_r = "
    f"{h_r:.4f} μm"
)
print(
    f"reflector 底部半徑 = "
    f"{r_profile[-1]:.4f} μm"
)
print()
print(
    "Step 6 條件是否全部滿足：",
    bool(np.all(condition_ok))
)
print(
    "Step 9 LED 包覆條件：",
    bool(step9_ok)
)
print("=" * 65)

plt.figure(
    figsize=(8, 6)
)

plt.plot(
    r_plot,
    Z_display,
    linewidth=2,
    label="Reflector profile"
)

plt.plot(
    -r_plot,
    Z_display,
    linewidth=2
)

plt.xlabel(
    "Radius R (μm)"
)

plt.ylabel(
    "Z (μm)"
)

plt.title(
    "2D Bowl Reflector Profile"
)

plt.grid(True)
plt.axis("equal")
plt.legend()
plt.tight_layout()

plt.figure(
    figsize=(8, 6)
)

plt.plot(
    z_profile + h_r,
    Sr_values,
    label="Sr design"
)

plt.plot(
    z_profile + h_r,
    Sr_min_values,
    "--",
    label="Sr lower bound"
)

plt.plot(
    z_profile + h_r,
    Sr_max_values,
    "--",
    label="Sr upper bound"
)

plt.xlabel(
    "Z (μm)"
)

plt.ylabel(
    "Sr = dr/dz"
)

plt.title(
    "Reflector Slope"
)

plt.grid(True)
plt.legend()
plt.tight_layout()

fig = plt.figure(
    figsize=(9, 7)
)

ax = fig.add_subplot(
    111,
    projection="3d"
)

ax.plot_surface(
    X3D,
    Y3D,
    Z3D,
    rstride=8,
    cstride=8,
    linewidth=0,
    antialiased=True
)

ax.set_xlabel(
    "X (μm)"
)

ax.set_ylabel(
    "Y (μm)"
)

ax.set_zlabel(
    "Z (μm)"
)

ax.set_title(
    "3D Bowl Reflector"
)

ax.set_box_aspect(
    (
        1,
        1,
        h_r
        /
        (2.0 * rr_design_0)
    )
)

plt.tight_layout()

plt.figure(
    figsize=(7, 7)
)

plt.scatter(
    X3D.flatten(),
    Y3D.flatten(),
    c=Z3D.flatten(),
    s=2
)

plt.xlabel(
    "X (μm)"
)

plt.ylabel(
    "Y (μm)"
)

plt.title(
    "Top View of Bowl Reflector"
)

plt.axis("equal")

plt.colorbar(
    label="Z (μm)"
)

plt.tight_layout()

if EXPORT_GRD:

    r_interp = r_plot.copy()
    z_interp = Z_display.copy()

    sort_index = np.argsort(
        r_interp
    )

    r_interp = (
        r_interp[
            sort_index
        ]
    )

    z_interp = (
        z_interp[
            sort_index
        ]
    )

    sag_function = interp1d(
        r_interp,
        z_interp,
        kind="linear",
        bounds_error=False,
        fill_value=(
            z_interp[0],
            z_interp[-1]
        )
    )

    R_top = r_interp[-1]

    grid_size = (
        2.4
        * R_top
    )

    x_grid = np.linspace(
        -grid_size / 2.0,
        grid_size / 2.0,
        grid_points
    )

    y_grid = np.linspace(
        -grid_size / 2.0,
        grid_size / 2.0,
        grid_points
    )

    X_grid, Y_grid = np.meshgrid(
        x_grid,
        y_grid
    )

    R_grid = np.sqrt(
        X_grid**2
        +
        Y_grid**2
    )

    Sag_um = sag_function(
        R_grid
    )

    R_bottom = (
        r_interp[0]
    )

    Sag_um[
        R_grid < R_bottom
    ] = z_interp[0]

    Sag_um[
        R_grid > R_top
    ] = z_interp[-1]

    Sag_mm = (
        Sag_um
        * 1e-3
    )

    dx_mm = (
        x_grid[1]
        - x_grid[0]
    ) * 1e-3

    dy_mm = (
        y_grid[1]
        - y_grid[0]
    ) * 1e-3

    try:
        script_dir = os.path.dirname(
            os.path.abspath(__file__)
        )
    except NameError:
        script_dir = os.getcwd()

    filename = os.path.join(
        script_dir,
        "Bowl_Reflector_From_Document.GRD"
    )

    with open(
        filename,
        "w"
    ) as f:

        f.write(
            f"{grid_points} "
            f"{grid_points} "
            f"{dx_mm:.9f} "
            f"{dy_mm:.9f} "
            f"0 0 0\n"
        )

        for j in range(
            grid_points
        ):

            for i in range(
                grid_points
            ):

                f.write(
                    f"{Sag_mm[j, i]:.9e} "
                    f"0 0 0\n"
                )

    print()
    print(
        "Zemax GRD 已輸出："
    )
    print(
        filename
    )

number_of_rays = 100000

reflectivity = 0.95

lambertian_exponent = 1.0

source_max_angle_deg = 90.0

source_max_angle = np.deg2rad(
    source_max_angle_deg
)

max_bounces = 5

rng = np.random.default_rng(
    seed=42
)

surface_r = interp1d(
    Z_display,
    r_plot,
    kind="linear",
    bounds_error=False,
    fill_value="extrapolate"
)

dr_dz_values = np.gradient(
    r_plot,
    Z_display
)

surface_dr_dz = interp1d(
    Z_display,
    dr_dz_values,
    kind="linear",
    bounds_error=False,
    fill_value="extrapolate"
)


def reflector_F(point):

    x_point = point[0]
    y_point = point[1]
    z_point = point[2]

    if (
        z_point < 0.0
        or
        z_point > h_r
    ):
        return np.nan

    rho = np.sqrt(
        x_point**2
        +
        y_point**2
    )

    reflector_radius = float(
        surface_r(z_point)
    )

    return (
        rho
        -
        reflector_radius
    )


def find_reflector_intersection(
    position,
    direction
):

    z0 = position[2]

    dz_ray = direction[2]

    if dz_ray > 0:

        t_limit = (
            (h_r - z0)
            /
            dz_ray
        )

    elif dz_ray < 0:

        t_limit = (
            (0.0 - z0)
            /
            dz_ray
        )

    else:

        t_limit = (
            4.0
            * max(
                ww_design,
                h_r
            )
        )

    if t_limit <= 1e-10:
        return None

    samples = 120

    t_values = np.linspace(
        1e-7,
        t_limit,
        samples
    )

    previous_t = t_values[0]

    previous_point = (
        position
        +
        previous_t
        * direction
    )

    previous_F = reflector_F(
        previous_point
    )

    for current_t in t_values[1:]:

        current_point = (
            position
            +
            current_t
            * direction
        )

        current_F = reflector_F(
            current_point
        )

        if (
            np.isnan(previous_F)
            or
            np.isnan(current_F)
        ):

            previous_t = current_t
            previous_F = current_F
            continue

        if (
            previous_F == 0
            or
            previous_F
            * current_F
            < 0
        ):

            def root_function(t):

                point = (
                    position
                    +
                    t
                    * direction
                )

                return reflector_F(
                    point
                )

            try:

                t_hit = brentq(
                    root_function,
                    previous_t,
                    current_t
                )

                hit_point = (
                    position
                    +
                    t_hit
                    * direction
                )

                return hit_point

            except ValueError:

                pass

        previous_t = current_t
        previous_F = current_F

    return None


def reflector_normal(point):

    x_point = point[0]
    y_point = point[1]
    z_point = point[2]

    rho = np.sqrt(
        x_point**2
        +
        y_point**2
    )

    if rho < 1e-12:

        return np.array(
            [
                0.0,
                0.0,
                1.0
            ]
        )

    slope = float(
        surface_dr_dz(
            z_point
        )
    )

    normal = np.array(
        [
            x_point / rho,
            y_point / rho,
            -slope
        ],
        dtype=float
    )

    normal = (
        normal
        /
        np.linalg.norm(
            normal
        )
    )

    return normal


def reflect_ray(
    direction,
    normal
):

    reflected_direction = (
        direction
        -
        2.0
        *
        np.dot(
            direction,
            normal
        )
        *
        normal
    )

    reflected_direction = (
        reflected_direction
        /
        np.linalg.norm(
            reflected_direction
        )
    )

    return reflected_direction


def check_waveguide(
    position,
    direction
):

    z0 = position[2]

    dz_ray = direction[2]

    if dz_ray <= 0:

        return (
            False,
            None,
            None
        )

    t_waveguide = (
        (h_r - z0)
        /
        dz_ray
    )

    if t_waveguide <= 0:

        return (
            False,
            None,
            None
        )

    hit_point = (
        position
        +
        t_waveguide
        * direction
    )

    x_hit = hit_point[0]
    y_hit = hit_point[1]

    rho_hit = np.sqrt(
        x_hit**2
        +
        y_hit**2
    )

    waveguide_radius = (
        ww_design
        /
        2.0
    )

    spatial_ok = (
        rho_hit
        <=
        waveguide_radius
    )

    cos_incident_angle = np.clip(
        direction[2]
        /
        np.linalg.norm(
            direction
        ),
        -1.0,
        1.0
    )

    incident_angle = np.arccos(
        cos_incident_angle
    )

    angular_ok = (
        incident_angle
        <=
        alpha
    )

    success = (
        spatial_ok
        and
        angular_ok
    )

    return (
        success,
        hit_point,
        incident_angle
    )


def generate_lambertian_ray():

    led_radius = (
        wl_design
        /
        2.0
    )

    u_position = rng.random()
    v_position = rng.random()

    source_radius = (
        led_radius
        *
        np.sqrt(
            u_position
        )
    )

    source_phi = (
        2.0
        *
        np.pi
        *
        v_position
    )

    x_source = (
        source_radius
        *
        np.cos(
            source_phi
        )
    )

    y_source = (
        source_radius
        *
        np.sin(
            source_phi
        )
    )

    z_source = 0.0

    position = np.array(
        [
            x_source,
            y_source,
            z_source
        ],
        dtype=float
    )

    u_angle = rng.random()
    v_angle = rng.random()

    sin_theta = (
        np.sin(
            source_max_angle
        )
        *
        np.sqrt(
            u_angle
        )
    )

    theta = np.arcsin(
        sin_theta
    )

    phi_angle = (
        2.0
        *
        np.pi
        *
        v_angle
    )

    direction = np.array(
        [
            np.sin(theta)
            *
            np.cos(phi_angle),

            np.sin(theta)
            *
            np.sin(phi_angle),

            np.cos(theta)
        ],
        dtype=float
    )

    direction = (
        direction
        /
        np.linalg.norm(
            direction
        )
    )

    return (
        position,
        direction,
        theta
    )


input_power = float(
    number_of_rays
)

direct_coupled_power = 0.0
reflected_coupled_power = 0.0
lost_power = 0.0

direct_ray_count = 0
reflected_ray_count = 0
lost_ray_count = 0

ray_paths = []

max_plot_rays = 200

accepted_angles = []

for ray_index in range(
    number_of_rays
):

    (
        position,
        direction,
        source_theta
    ) = generate_lambertian_ray()

    ray_power = 1.0

    bounce_count = 0

    path = [
        position.copy()
    ]

    coupled = False

    while (
        bounce_count
        <=
        max_bounces
    ):

        reflector_hit = (
            find_reflector_intersection(
                position,
                direction
            )
        )

        (
            waveguide_success,
            waveguide_hit,
            incident_angle
        ) = check_waveguide(
            position,
            direction
        )

        t_waveguide = np.inf

        if (
            direction[2] > 0
            and
            waveguide_hit is not None
        ):

            displacement = (
                waveguide_hit
                -
                position
            )

            t_waveguide = np.linalg.norm(
                displacement
            )

        t_reflector = np.inf

        if reflector_hit is not None:

            displacement = (
                reflector_hit
                -
                position
            )

            t_reflector = np.linalg.norm(
                displacement
            )

        if (
            t_waveguide
            <
            t_reflector
        ):

            path.append(
                waveguide_hit.copy()
            )

            if waveguide_success:

                coupled = True

                accepted_angles.append(
                    np.rad2deg(
                        incident_angle
                    )
                )

                if bounce_count == 0:

                    direct_coupled_power += (
                        ray_power
                    )

                    direct_ray_count += 1

                else:

                    reflected_coupled_power += (
                        ray_power
                    )

                    reflected_ray_count += 1

                if (
                    len(ray_paths)
                    <
                    max_plot_rays
                ):

                    ray_paths.append(
                        np.array(
                            path
                        )
                    )

            else:

                lost_power += (
                    ray_power
                )

                lost_ray_count += 1

            break

        elif reflector_hit is not None:

            path.append(
                reflector_hit.copy()
            )

            normal = (
                reflector_normal(
                    reflector_hit
                )
            )

            direction = (
                reflect_ray(
                    direction,
                    normal
                )
            )

            ray_power *= (
                reflectivity
            )

            bounce_count += 1

            position = (
                reflector_hit
                +
                direction
                *
                1e-6
            )

            if ray_power < 1e-8:

                lost_power += (
                    ray_power
                )

                lost_ray_count += 1

                break

        else:

            lost_power += (
                ray_power
            )

            lost_ray_count += 1

            break

    if (
        not coupled
        and
        bounce_count
        >
        max_bounces
    ):

        lost_power += (
            ray_power
        )

        lost_ray_count += 1

direct_efficiency = (
    direct_coupled_power
    /
    input_power
)

reflected_efficiency = (
    reflected_coupled_power
    /
    input_power
)

total_coupling_efficiency = (
    direct_efficiency
    +
    reflected_efficiency
)

coupling_loss_db = (
    -10.0
    *
    np.log10(
        max(
            total_coupling_efficiency,
            1e-15
        )
    )
)

print()
print("=" * 70)
print(
    "Lambertian Source Coupling Efficiency Test"
)
print("=" * 70)

print(
    f"Lambertian cosine exponent = "
    f"{lambertian_exponent:.1f}"
)

print(
    f"Source angular range = "
    f"0 ~ {source_max_angle_deg:.1f} deg"
)

print(
    f"Number of rays = "
    f"{number_of_rays}"
)

print(
    f"Reflector reflectivity = "
    f"{reflectivity:.4f}"
)

print(
    f"Waveguide acceptance angle = "
    f"{np.rad2deg(alpha):.4f} deg"
)

print("-" * 70)

print(
    f"Direct coupled rays = "
    f"{direct_ray_count}"
)

print(
    f"Reflected coupled rays = "
    f"{reflected_ray_count}"
)

print(
    f"Lost rays = "
    f"{lost_ray_count}"
)

print("-" * 70)

print(
    f"Direct coupling efficiency = "
    f"{direct_efficiency * 100:.4f} %"
)

print(
    f"Reflector coupling efficiency = "
    f"{reflected_efficiency * 100:.4f} %"
)

print(
    f"TOTAL coupling efficiency = "
    f"{total_coupling_efficiency * 100:.4f} %"
)

print(
    f"Coupling loss = "
    f"{coupling_loss_db:.4f} dB"
)

print("=" * 70)

plt.figure(
    figsize=(9, 7)
)

plt.plot(
    r_plot,
    Z_display,
    linewidth=2,
    label="Reflector"
)

plt.plot(
    -r_plot,
    Z_display,
    linewidth=2
)

plt.plot(
    [
        -wl_design / 2.0,
        wl_design / 2.0
    ],
    [
        0.0,
        0.0
    ],
    linewidth=4,
    label="Lambertian LED"
)

plt.plot(
    [
        -ww_design / 2.0,
        ww_design / 2.0
    ],
    [
        h_r,
        h_r
    ],
    linewidth=4,
    label="Waveguide"
)

for ray_path in ray_paths:

    plt.plot(
        ray_path[:, 0],
        ray_path[:, 2],
        linewidth=0.7,
        alpha=0.4
    )

plt.xlabel(
    "X (μm)"
)

plt.ylabel(
    "Z (μm)"
)

plt.title(
    "Lambertian Ray Tracing"
)

plt.grid(True)
plt.legend()
plt.axis("equal")
plt.tight_layout()

if len(
    accepted_angles
) > 0:

    plt.figure(
        figsize=(8, 6)
    )

    plt.hist(
        accepted_angles,
        bins=40
    )

    plt.axvline(
        np.rad2deg(alpha),
        linestyle="--",
        label="Waveguide acceptance angle"
    )

    plt.xlabel(
        "Incident Angle (deg)"
    )

    plt.ylabel(
        "Ray Count"
    )

    plt.title(
        "Accepted Ray Incident Angle Distribution"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

plt.show()