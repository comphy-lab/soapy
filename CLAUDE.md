# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# SOAPY - Simulation Of Axisymmetric Phenomena in bubblY films

High-fidelity CFD simulations of soap bubble dynamics studying Taylor-Culick instabilities, wrinkling patterns, and film drainage using the Basilisk C framework.

## Architecture Overview

SOAPY extends the Basilisk CFD library to simulate two-phase flows in thin liquid films. The codebase consists of:

1. **Basilisk Core** (`basilisk/src/`): The underlying CFD framework providing:
   - Adaptive mesh refinement (AMR) infrastructure
   - Volume-of-Fluid (VOF) interface tracking
   - Navier-Stokes solvers with surface tension
   - Axisymmetric coordinate transformations

2. **Simulation Cases** (`simulationCases/`): Self-contained simulation programs that:
   - Model complete bubble geometry (spherical cap + toroidal rim)
   - Track smoke tracer diffusion for flow visualization
   - Implement adaptive refinement near interfaces
   - Output binary snapshots for post-processing

3. **Post-Processing Tools** (`postProcess/`):
   - `getData.c`: Extracts field data from binary dumps
   - `getFacet.c`: Reconstructs interface geometry from VOF field
   - `Video-generic.py`: Generates animations and analysis plots

## Build System

### Initial Setup (One-time)
```bash
# Install Basilisk locally and configure environment
./reset_install_requirements.sh
source .project_config  # Sets BASILISK path
```

### Compilation Commands
```bash
# Standard simulation compilation
cd simulationCases
make soapBubble-full.tst   # Creates executable with .tst extension

# Manual compilation with qcc (Basilisk's compiler wrapper)
qcc -O2 -disable-dimensions -autolink soapBubble-full.c -o soapBubble-full -lm

# Post-processing tools
cd postProcess
qcc -autolink getData.c -o getData -lm
qcc -autolink getFacet.c -o getFacet -lm
```

### Running Simulations
```bash
# Execute simulation (outputs to intermediate/ directory)
./soapBubble-full.tst

# Key output files:
# - dump-* : Binary snapshots at regular intervals
# - log : Time series data (energies, extrema)
# - out : Standard output/debug information
```

## Key Basilisk Headers Used

The project primarily uses these Basilisk modules:

- **`axi.h`**: Enables axisymmetric (2D cylindrical) coordinates
- **`navier-stokes/centered.h`**: Incompressible Navier-Stokes solver
- **`navier-stokes/conserving.h`**: Conservative momentum advection for VOF
- **`two-phase.h`**: Two-phase flow with density/viscosity ratios
- **`tension.h`**: Surface tension using Continuous Surface Force (CSF)
- **`henry.h`**: Gas dissolution/Henry's law implementation
- **`tracer.h`**: Passive scalar advection (smoke visualization)
- **`diffusion.h`**: Scalar diffusion solver

## Simulation Structure Pattern

All simulation files follow this structure:

1. **Documentation header**: Comprehensive markdown documentation
2. **Include directives**: Basilisk headers in specific order
3. **Macro definitions**: Physical parameters, tolerances
4. **Boundary conditions**: Using Basilisk's boundary syntax
5. **Global variables**: Simulation parameters, output controls
6. **Main function**: Parameter initialization, `run()` call
7. **Event functions**:
   - `event init()`: Initial conditions
   - `event adapt()`: Adaptive refinement criteria
   - `event writingFiles()`: Data output
   - `event logWriting()`: Time series logging

## Code Style Guidelines

- **Indentation**: 2 spaces (no tabs)
- **Line length**: 80 characters maximum
- **Documentation**: Markdown in comments starting with `/**`
- **Variable naming**: snake_case (e.g., `max_level`)
- **Function naming**: camelCase (e.g., `writingFiles`)
- **Physical parameters**: UPPERCASE macros (e.g., `Oh1`, `PECLETGAS`)
- **Error tolerances**: Scientific notation (e.g., `1e-6`)

## Common Development Tasks

### Modifying Physical Parameters
Edit macro definitions in the simulation file:
```c
#define Oh1 (1e-3)        // Ohnesorge number
#define k (25)            // R/h ratio
#define MAXlevel (10)     // Maximum refinement level
```

### Adding New Diagnostics
Add to the `logWriting` event:
```c
event logWriting (t += 0.01; t <= tmax) {
  scalar speed[];
  foreach()
    speed[] = norm(u);
  stats s = statsf(speed);
  fprintf(ferr, "%g %g %g %g\n", t, s.min, s.max, s.sum);
}
```

### Extracting Data from Dumps
```bash
# Extract velocity field at t=0.5
./getData dump-0.5 -v u.x u.y > velocity_t0.5.dat

# Extract interface facets
./getFacet dump-0.5 > interface_t0.5.dat
```

## Important Notes

- The `src-local/` directory mentioned in compile commands does not currently exist
- Basilisk uses `qcc` which wraps the system C compiler (clang on macOS)
- The `-disable-dimensions` flag turns off dimensional analysis checks
- Binary dumps can be large; the `intermediate/` directory needs adequate space
- Post-processing typically requires Python with matplotlib, numpy installed