# Bowl Reflector Design and Ray-Tracing Simulation

A Python research prototype for designing an axisymmetric bowl-shaped reflector, generating a sampled surface for optical simulation, and estimating geometric light collection with ray tracing.

>  ## Project Status

The optical component structure design and Python-based numerical validation have been completed, achieving a simulated coupling efficiency of **95.6774%**.

Further improvements to the bowl-shaped reflector geometry, particularly extending the reflector to fully surround the light source, are expected to potentially increase the coupling efficiency to **over 98%**. However, this improvement remains a design hypothesis and has not yet been verified.

In **Ansys Zemax OpticStudio**, the optical system achieved a simulated coupling efficiency of approximately **90%** with a 30° collection half-angle. However, reducing the collection half-angle to 15° resulted in a decrease in coupling efficiency to approximately **79%**.

The discrepancy between the Python-based simulation and the Zemax results, as well as the sensitivity of coupling efficiency to the collection half-angle, is currently under investigation. 

## Features

- Reflector profile calculated by numerical integration (`scipy.integrate.solve_ivp`)
- Geometric constraint checks and 2D/3D visualization
- Export of a sampled surface to `Bowl_Reflector_From_Document.GRD` for subsequent Zemax-format evaluation 
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


## Key Model Settings

The following parameters are defined in `main.py` and used for reflector geometry design and Monte Carlo ray-tracing simulation.

| Parameter | Value | Description |
|---|---|---|
| LED design width (`wl_design`) | 10 μm | Diameter of the modeled circular light source |
| Waveguide design width (`ww_design`) | 50 μm | Diameter of the modeled circular waveguide acceptance region |
| Waveguide acceptance half-angle (`alpha_deg`) | 30° | Maximum accepted ray angle relative to the optical axis |
| Reflector depth ratio | 0.90 | Ratio used to determine the reflector depth |
| Reflectivity | 0.95 | Assumed power reflectivity per reflection |
| Maximum reflections | 5 | Configured reflection limit |
| Number of rays | 100,000 | Monte Carlo simulation sample size |
| Random seed | 42 | Reproducible random sampling |
| Reflector profile points | 800 | Numerical profile sampling |
| Azimuthal points | 240 | 3D reflector surface sampling |
| GRD grid resolution | 201 × 201 | Sampled surface resolution for Zemax |

The current ray-direction sampling corresponds to a Lambertian cosine-weighted distribution over the forward hemisphere.

## Model limitations

- The ray-tracing acceptance rule checks position and angle at a target plane; it does not simulate full waveguide modes or wave-optical coupling.
- The `lambertian_exponent` variable is printed but is not used to parameterize the ray-direction sampling in the current implementation.
- The generated `.GRD` file has been successfully imported into Ansys Zemax OpticStudio for optical simulation. However, the current simulation results have not yet met the expected performance targets, and further investigation and optimization are ongoing.
- The power-accounting and intersection logic should be validated before quantitative results are treated as reliable.


## Results

### 1. Bowl Reflector Geometry

The axisymmetric bowl-shaped reflector was generated using numerical integration of the reflector profile equation.

| Design Parameter | Result |
|---|---|
| Initial reflector width | 47.0000 μm |
| Initial reflector radius | 23.5000 μm |
| Calculated design height | 53.9615 μm |
| Actual reflector depth | 48.5654 μm |
| Bottom reflector radius | 6.2309 μm |
| Step 6 geometric condition | Passed |
| Step 9 LED enclosure condition | Not satisfied |

The reflector geometry satisfies the implemented Step 6 geometric condition. However, the current design does not fully satisfy the LED enclosure condition.

Further geometric optimization is planned to improve light collection by extending the reflector around the source.

### 2. Python Monte Carlo Ray-Tracing Results

A Monte Carlo ray-tracing simulation was performed using 100,000 rays emitted from a circular Lambertian source.

The model includes specular reflection, an assumed reflectivity of 95% per bounce, and position-based and angle-based acceptance checks at the waveguide plane.

| Simulation Metric | Result |
|---|---|
| Number of simulated rays | 100,000 |
| Directly accepted rays | 18,745 |
| Reflected and accepted rays | 81,240 |
| Rays classified as lost | 15 |
| Direct coupling efficiency | 18.7450% |
| Reflected coupling efficiency | 76.9324% |
| **Total simulated coupling efficiency** | **95.6774%** |
| **Estimated coupling loss** | **0.1919 dB** |

The reported efficiencies are based on ray-power accounting, including the modeled reflection losses. The total coupling efficiency is calculated as the sum of the directly accepted and reflected accepted power, divided by the initial emitted power.

These results represent geometric ray-tracing performance under the current model assumptions rather than experimentally measured or waveguide-mode-overlap coupling efficiency.

### 3. Ansys Zemax OpticStudio Simulation

The reflector surface was exported as `Bowl_Reflector_From_Document.GRD` and successfully imported into Ansys Zemax OpticStudio for optical simulation.

The current reported Zemax results are:

| Collection Half-Angle | Reported Efficiency |
|---|---|
| 30° | Approximately 90% |
| 15° | Approximately 79% |

The Zemax results are lower than the Python-based simulation result of 95.6774% under the reported 30° acceptance condition.

The differences between the Python and Zemax simulation results are currently under investigation, including the effects of collection angle, optical geometry, and simulation settings.

### 4. Future Optimization

The current reflector design has not fully satisfied the LED enclosure condition.

Further optimization will investigate whether extending the reflector geometry to fully surround the light source can improve simulated coupling efficiency beyond 98%.

This remains a proposed design target and has not yet been demonstrated.

Additional plots, ray-tracing visualizations, and simulation screenshots will be organized in the `results/` directory.

Existing Zemax simulation screenshots and measurements are available in [`docs/Zemax_Simulation_Results.docx`](docs/Zemax_Simulation_Results.docx).


## Author

**Jian-Yang Zhao**  
M.S. Student, Department of Electronic and Computer Engineering  
National Taiwan University of Science and Technology (Taiwan Tech)

**Research Interests:** Optical Interconnects, Optical System Design, Ray-Tracing Simulation, and Python-Based Numerical Modeling.

**GitHub:** [ZHAOJIAN-YANG](https://github.com/ZHAOJIAN-YANG)


## Existing Zemax results

The current Zemax simulation screenshots and measurements are archived in [`docs/Zemax_Simulation_Results.docx`](docs/Zemax_Simulation_Results.docx). Additional result figures and validation notes will be added later. The current Zemax simulation evaluates optical power collection at the detector plane. Further analysis of waveguide mode matching and coupling characteristics will be conducted to assess the overall optical coupling performance.
