/**
# Die Injection for Flow Visualization

This module introduces a circular tracer (die) into the flow at a specified time 
and location. The die then advects with the flow, allowing visualization of flow 
patterns. This is useful for visualizing complex flow structures in simulations
like lid-driven cavity flow.

## Parameters
- `tInjection`: time at which to inject the die
- `xInjection`, `yInjection`: position where the die is injected
- `dieRadius`: radius of the circular die
*/

#include "tracer.h"

// Die tracer parameters (can be overridden by the user)
double tInjection = 0.1;  // Default injection time
double xInjection = 0.0;  // Default X-position for injection
double yInjection = 0.0;  // Default Y-position for injection
double dieRadius = 0.05;  // Default radius of the circular die

// Define the scalar field for the die
scalar T[];
scalar * tracers = {T};

// Initialize the die tracer to zero everywhere
event init (t = 0) {
  foreach()
    T[] = 0.0;
}

// Inject the die at the specified time
event inject_die (t = tInjection) {
  fprintf(stderr, "Injecting die at t = %g, position = (%g, %g), radius = %g\n", 
          t, xInjection, yInjection, dieRadius);
  
  // Set die concentration to 1.0 within the circular region
  foreach() {
    double dist = sqrt(sq(x - xInjection) + sq(y - yInjection));
    if (dist <= dieRadius)
      T[] = 1.0;
  }
}