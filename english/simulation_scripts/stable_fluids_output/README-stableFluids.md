# Stable Fluids Implementation Comparison: Python vs Fortran

## Key Differences Identified

### 1. Grid Type Fundamental Difference
- **Python**: Node-centered (vertex-based) grid
- **Fortran**: Cell-centered (staggered) grid

This is the most fundamental difference that affects all calculations.

### 2. Grid Resolution
- **Python**: 41x41 nodes (N_POINTS = 41)
- **Fortran**: 42x42 cells with ghost cells (igx = jgx = 42)

### 3. Grid Spacing
- **Python**: `element_length = 1.0 / (41-1) = 0.025`
- **Fortran**: `dxg0 = 1.0 / (42-2) = 0.025`

While the spacing appears the same, they represent different things (node spacing vs cell width).

### 4. Grid Origin and Boundaries
- **Python**: 
  - Grid starts exactly at (0.0, 0.0)
  - Uses `np.linspace(0.0, 1.0, 41)` for uniform node distribution
  - Boundaries are at exact domain edges

- **Fortran**: 
  - Grid origin: `xmg(0) = -dxg0*0.5 = -0.0125`
  - First cell center is at approximately (-0.0125, -0.0125)
  - Uses ghost cells for boundary conditions

### 5. Boundary Condition Implementation
- **Python**: 
  - Implicit through differential operators
  - Uses array slicing `[1:-1, 1:-1]` to exclude boundaries
  - Derivatives automatically set to zero at boundaries

- **Fortran**: 
  - Explicit boundary conditions with ghost cells (margin=1)
  - Separate boundary condition subroutines
  - Ghost cells updated at each step

### 6. Differential Operators
- **Python**: 
  ```python
  diff[1:-1, 1:-1] = (field[2:, 1:-1] - field[0:-2, 1:-1]) / (2 * element_length)
  ```
  - Central differences on interior nodes only
  - Boundary values remain zero

- **Fortran**: 
  - Operates on all cells including ghost cells
  - Different indexing due to cell-centered approach

### 7. Advection Implementation
- **Python**: 
  - Uses scipy's `interpolate.interpn` for backtracing
  - Interpolates on node values directly

- **Fortran**: 
  - Custom interpolation routine
  - Must handle cell-centered values differently

### 8. Linear Solvers
- **Python**: 
  - Conjugate Gradient (CG) solver from scipy
  - No iteration limit (MAX_ITER_CG = None)
  - Matrix-free implementation using LinearOperator
  - Generally faster convergence

- **Fortran**: 
  - Jacobi iteration with MAX_ITER_CG = 100
  - Tolerance-based convergence (1e-6)
  - Updates solution in-place
  - Slower convergence but simpler implementation

## Impact on Results

These differences lead to:

1. **Different numerical accuracy**: Cell-centered vs node-centered grids have different truncation errors
2. **Different boundary behavior**: Ghost cells vs implicit boundaries affect flow near walls
3. **Different interpolation**: Advection step will produce different results
4. **Different convergence**: Jacobi vs CG have different convergence characteristics

## Recommendations for Matching Results

To make the implementations more comparable:

1. **Standardize grid type**: Either convert Fortran to node-centered or Python to cell-centered
2. **Match grid resolution**: Use same number of computational points
3. **Align boundary handling**: Implement consistent boundary conditions
4. **Use same solver**: Implement CG in Fortran or Jacobi in Python
5. **Verify interpolation**: Ensure advection schemes are equivalent

## Visualization Differences

The different grid types also affect visualization:
- Python plots show values at nodes
- Fortran plots show cell-averaged values
- This can make direct visual comparison misleading

## Specific Implementation Differences

### Laplacian Operator
**Python:**
```python
diff[1:-1, 1:-1] = (field[0:-2, 1:-1] + field[1:-1, 0:-2] - 4*field[1:-1, 1:-1] 
                    + field[2:, 1:-1] + field[1:-1, 2:]) / element_length**2
```

**Fortran:** 
- Uses separate `laplacian` subroutine with ghost cell handling
- Different indexing due to cell-centered approach

### Divergence Calculation
**Python:**
```python
divergence = partial_derivative_x(u) + partial_derivative_y(v)
```

**Fortran:**
- Explicit divergence calculation with proper ghost cell treatment
- Central differences computed on cell edges

### Time Step Differences
- **Python**: TIME_STEP_LENGTH = 0.1 (fixed)
- **Fortran**: dt0 = 0.1 (can be adaptive based on CFL condition)

## Conclusion

The implementations are solving the same physical equations but with fundamentally different numerical discretizations. Perfect agreement is not expected without aligning these core differences. The main sources of discrepancy are:

1. **Grid discretization**: Node-centered vs cell-centered
2. **Boundary treatment**: Implicit vs explicit ghost cells
3. **Numerical methods**: Different solver implementations
4. **Grid resolution**: 41x41 nodes vs 42x42 cells

To achieve comparable results, one would need to either:
- Convert both to the same grid type
- Implement identical boundary conditions
- Use the same linear solvers
- Ensure grid resolutions match properly