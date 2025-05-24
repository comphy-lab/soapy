# Basilisk-101 Development Guidelines

## project structure 
├-- `basilisk/src/`: Core Basilisk CFD library (reference only, do not modify) <- This contains the codebase for Basilisk C--a wrapper around C to do computational fluid dynamics. 
## Some selected files inside basilisk/src:
│   ├── `navier-stokes/centered.h`: Main centered Navier-Stokes solver
│   ├── `navier-stokes/conserving.h`: Conservative form solver with VoF momentum advection
│   ├── `two-phase*.h`: Two-phase flow implementations (VoF/Level-set/CLSVOF)
│   ├── `reduced.h`: this implements the reduced gravity approach. 
│   ├── `curvature.h`: this implements the curvature module for computing interface properties.
│   ├── `tension.h`: this implements the surface tension force using the Brackbill method.
│   ├── `integral.h`: Integral formulation for surface tension.
│   ├── `tracer.h`: this implements the advection equation for passive tracers.
│   ├── `diffusion.h`: this implements the diffusion equation for passive tracers.
│   ├── `vof.h`: this implements the volume of fluid (VOF) method.
│   ├── `viscosity.h`: this implements the implicit viscous stress solver in Basilisk.
│   ├── `axi.h`: this implements axisymmetric metric so that 2D equations can be extended to axi. 
├── `postProcess/`: Project-specific post-processing tools
├── `src-local/`: Custom header files extending Basilisk functionality
├── `testCases/`: Test cases with their own Makefile

## Code Style
- Indentation: 2 spaces (no tabs)
- Line length: 80 characters max
- Documentation: Use markdown in comments starting with `/**`. Do not use `*` in comments.
- Spacing: Space after commas, spaces around operators (+, -)
- Files: Core functionality in `.h` headers, tests in `.c` files
- Naming: Snake_case for variables, camelCase for functions
- Error handling: Return values with stderr messages

## Build & Test Commands
- Compile: `qcc -autolink file.c -o executable -lm`
- Compile with specific headers: `qcc -I$PWD/src-local -autolink file.c -o executable -lm`