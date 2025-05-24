# SOAPY - Simulation Of Axisymmetric Phenomena in bubblY films

[![Issues](https://img.shields.io/github/issues/comphy-lab/SOAPY)](https://github.com/comphy-lab/SOAPY/issues)
[![License](https://img.shields.io/github/license/comphy-lab/SOAPY)](https://github.com/comphy-lab/SOAPY/blob/main/LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/comphy-lab/SOAPY)](https://github.com/comphy-lab/SOAPY/commits/main)
[![Basilisk](https://img.shields.io/badge/Basilisk-Compatible-green)](http://basilisk.fr/)
[![Research](https://img.shields.io/badge/Research-Fluid%20Dynamics-blue)](https://comphy-lab.org)
![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/comphy-lab/SOAPY?utm_source=oss&utm_medium=github&utm_campaign=comphy-lab%2FSOAPY&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)

**S**imulation **O**f **A**xisymmetric **P**henomena in bubbl**Y** films

High-fidelity simulations of soap bubble dynamics and wrinkling instabilities using Basilisk C.

## Research Overview

This repository contains computational tools for studying the dynamics of thin liquid films forming bubbles, with a focus on wrinkling instabilities driven by surface tension and viscous effects. The simulations solve the two-phase Navier-Stokes equations in axisymmetric configurations, capturing the complex interplay between:

- **Surface tension-driven dynamics** at the liquid-gas interface
- **Viscous dissipation** within the thin film
- **Gravitational effects** characterized by the Bond number
- **Film thickness variations** leading to instability formation

### Physical System

The code models a soap bubble consisting of:
- A thin liquid film (water with surfactants)
- Surrounding gas phase (air)
- Interface dynamics tracked using Volume-of-Fluid (VOF) method
- Smoke tracer diffusion for flow visualization

### Key Dimensionless Parameters

- **Ohnesorge number (Oh)**: Ratio of viscous to inertio-capillary forces
  $$\text{Oh} = \frac{\mu}{\sqrt{\rho\sigma R}}$$
- **Bond number (Bo)**: Ratio of gravitational to surface tension forces
  $$\text{Bo} = \frac{\rho g R^2}{\sigma}$$
- **Péclet number (Pe)**: Ratio of advective to diffusive transport
  $$\text{Pe} = \frac{V_\gamma R}{D}$$
  where $V_\gamma$ is the inertio-capillary velocity
- **Curvature parameter (k)**: Ratio of bubble radius to film thickness
  $$k = \frac{R}{h}$$

## Repository Structure

```
├── basilisk/src/               Core Basilisk CFD library (reference only, do not modify)
├── simulationCases/                Test cases for simulation
│   └── soapBubble-full.c           Full bubble simulation (spherical cap + rim)
│   └── soapBubble-half.c           Half-bubble simulation variant
└── postProcess/                Project-specific post-processing tools
    ├── getData.c            Data extraction utility
    ├── getFacet.c            Interface facet extraction
    └── Video-generic.py     Python script for post-processing
```

## Simulation Configurations

### Full Bubble Simulation (`soapBubble-full.c`)
Models a complete bubble with:
- Lower hemisphere: Circular cap geometry
- Upper region: Spherical rim structure
- Initial film thickness: $h = 1/k$
- Adaptive mesh refinement near interfaces

### Half Bubble Simulation (`soapBubble-half.c`)
Variant configuration for studying partial bubble dynamics with modified boundary conditions.

## Physical Parameters Used

| Parameter | Symbol | Typical Value | Description |
|-----------|--------|---------------|-------------|
| Ohnesorge number | $\text{Oh}$ | $10^{-3}$ | Viscous effects |
| Density ratio | $\rho_2/\rho_1$ | $10^{-3}$ | Gas/liquid density |
| Viscosity ratio | $\mu_2/\mu_1$ | $10^{-3}$ | Gas/liquid viscosity |
| Curvature parameter | $k$ | 25 | $R/h$ ratio |
| Péclet number (gas) | $\text{Pe}_{\text{gas}}$ | $10^{-1}$ | Smoke diffusion |
| Domain size | $L$ | 5.0 | Computational domain |

## Installation and Usage

### Prerequisites
- C99-compliant compiler (gcc or clang)
- GNU make
- Python 3+ with matplotlib, numpy for post-processing

### Quick Setup
```bash
# Clone repository
git clone git@github.com:comphy-lab/SOAPY.git
cd SOAPY

# Install Basilisk locally
./reset_install_requirements.sh
source .project_config

# Compile simulations
cd simulationCases
make soapBubble-full.tst
```

### Running Simulations
```bash
# Run full bubble simulation
./soapBubble-full.tst

# Run with custom parameters (edit source file for parameter changes)
# Key parameters: Oh1, k, tmax, MAXlevel

# Output files are saved in intermediate/ directory
```

### Post-Processing
```bash
cd postProcess

# Compile data extraction tools
qcc -autolink getData.c -o getData -lm
qcc -autolink getFacet.c -o getFacet -lm

# Generate animations
python Video-generic.py --case ../simulationCases --folder output_frames
```

## Output Data

The simulations generate:
- **Snapshot files**: Binary dumps at regular intervals ($\Delta t = 0.01$)
- **Interface data**: VOF field and curvature information
- **Velocity fields**: Full velocity components in axisymmetric coordinates
- **Tracer concentration**: Smoke diffusion patterns
- **Pressure fields**: Dynamic pressure evolution

## Scientific Background

The simulations investigate the Taylor-Culick instability in soap bubbles, where:

- **Thin film dynamics**: The bubble consists of a liquid film of thickness $h \ll R$ (bubble radius)
- **Wrinkling mechanism**: Instabilities arise from the competition between surface tension and viscous forces
- **Rim formation**: The bubble edge develops a toroidal rim structure that influences stability
- **Drainage patterns**: Gravity and capillary forces create complex flow patterns within the film

### Relevant Physics

- **Taylor-Culick velocity**: Characteristic retraction speed
  $$V_{\text{TC}} = \sqrt{\frac{2\sigma}{\rho h}}$$
- **Rayleigh-Plateau instability**: Drives rim breakup and droplet formation
- **Marangoni effects**: Surface tension gradients due to surfactant concentration
- **Film thinning**: Coupled drainage and evaporation processes

## Numerical Methods

The code employs:

- **Volume-of-Fluid (VOF)** method for interface tracking
- **Adaptive mesh refinement** (AMR) for resolving thin films
- **Implicit viscous solver** for stability at low Oh numbers
- **Conservative momentum advection** for accurate interface dynamics
- **Height-function curvature** computation for surface tension

## Contributing

This code is part of ongoing research. For collaborations or questions:

- Open an issue for bug reports or feature requests
- Contact the authors for research collaborations
- See CLAUDE.md for code style guidelines

## Authors

- **Saumili Jana** - Primary developer (<jsaumili@gmail.com>)
- **Vatsal Sanjay** - [vatsalsy@comphy-lab.org](mailto:vatsalsy@comphy-lab.org)

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.
