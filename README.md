# Bowl Reflector Design and Ray-Tracing Simulation

A Python research prototype for designing an axisymmetric bowl-shaped reflector, generating a sampled surface for optical simulation, and estimating geometric light collection with ray tracing.

> **Status:** Research prototype. The simulation results are model-dependent and are not experimental validation. Review permission to share research code and data before making this repository public.

## Features

- Reflector profile calculated by numerical integration (`scipy.integrate.solve_ivp`)
- Geometric constraint checks and 2D/3D visualization
- Export of a sampled surface to `Bowl_Reflector_From_Document.GRD` for subsequent Zemax-format evaluation (import compatibility not independently verified)
- Monte Carlo ray tracing from a circular source, specular reflections, and geometric acceptance checks at the waveguide plane
- Reported direct/reflected coupling efficiency and estimated coupling loss

## Requirements

- Python 3.10+ (suggested; not independently tested across versions)
- NumPy, SciPy, Matplotlib

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

## Run

```bash
python main.py
```

The script displays plots, prints diagnostic values, and (when `EXPORT_GRD = True`) writes a `.GRD` file beside `main.py`.

**Runtime note:** The default setting traces 100,000 rays. For a quick first run, temporarily change `number_of_rays = 100000` to `number_of_rays = 1000`.

## Key model settings

The settings are defined near the beginning of `main.py` and in the ray-tracing section:

- Design widths and tolerances in µm
- `alpha_deg = 30.0`: modeled waveguide angular acceptance threshold
- `reflector_depth_ratio = 0.90`
- `reflectivity = 0.95`: assumed reflectivity per bounce
- `max_bounces = 5`
- `number_of_rays = 100000`
- Random generator seed: 42

## Model limitations

- The ray-tracing acceptance rule checks position and angle at a target plane; it does not simulate full waveguide modes or wave-optical coupling.
- The `lambertian_exponent` variable is printed but is not used to parameterize the ray-direction sampling in the current implementation.
- Exporting a `.GRD` file does not by itself verify correct import or optical performance in Zemax.
- The power-accounting and intersection logic should be validated before quantitative results are treated as reliable.

## Results

Add approved screenshots of the 2D profile, 3D reflector, ray tracing, and numerical output to `results/`. Include the parameter settings used for each figure.

## Author / attribution

Add your name, affiliations if authorized, and relevant citations before publishing.


## Existing Zemax results

The current Zemax simulation screenshots and measurements are archived in [`docs/Zemax_Simulation_Results.docx`](docs/Zemax_Simulation_Results.docx). Additional result figures and validation notes will be added later. Detector power readings are not automatically equivalent to guided-mode coupling efficiency.
