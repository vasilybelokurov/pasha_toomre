# Simulation Library Implementation Plan

1. Create the package metadata and immutable simulation configuration.
2. Implement the documented fourth-order Plummer acceleration and its registry.
3. Reconstruct the equations of motion and replace the vector-valued crossing event with a scalar terminal event.
4. Implement the two-pass DOP853 workflow, dense-output evaluation, and adaptive report grid.
5. Implement the 12 common diagnostics plus the extended intruder-relative `eta` diagnostic.
6. Write text, CSV, extrema, and JSON outputs.
7. Reconstruct the two-figure overview layout and the four-figure extended layout.
8. Add a shared CLI and two thin Venus driver scripts.
9. Add unit and smoke tests for forces, events, integration, diagnostics, and conservation.
10. Run and debug the test suite with the project scientific Python environment.
11. Execute both complete Venus runs through `tau=80`.
12. Inspect all reports and plots, then document verified commands and numerical results.

The original Word documents and JPEG remain unchanged.
