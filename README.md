# ternary-electromagnetism

EM field simulation on ternary lattices — electric and magnetic fields, Maxwell's equations, and wave propagation where field values live in **{-1, 0, +1}**.

---

## What It Is

Classical electromagnetism operates over continuous real-valued fields. This crate discretizes those fields onto a **ternary lattice** — each field value is restricted to one of three states: negative (-1), zero (0), or positive (+1). This discretization is natural for:

- **Digital EM hardware** where three-state logic gates drive field computations
- **Trit-based computing** architectures where balanced ternary arithmetic is native
- **Topological field theories** where charge quantization forbids intermediate values
- **Constraint-based simulation** where field states encode logical relationships rather than magnitudes

The library implements the full set of classical EM primitives — Coulomb's law, Biot-Savart, Maxwell's curl equations, wave propagation, polarization, and interference — all operating in the ternary regime.

---

## The Yee Lattice and Ternary Fields

The **Yee lattice** (introduced by Kane Yee in 1966) is a staggered-grid discretization of Maxwell's equations in which E and B fields are offset by half a cell in both space and time. This staggering gives the scheme several remarkable properties:

- **Second-order accuracy** in both space and time
- **Exact discrete charge conservation** — the discrete divergence of E is preserved
- **No spurious modes** — the leapfrog time-stepping is symplectic and energy-conserving

In the ternary setting, the Yee lattice is especially natural. The curl operator on a ternary lattice produces values in {-2, -1, 0, +1, +2}, which clamp back to {-1, 0, +1} under the update rule. This makes the ternary Yee scheme a **closed system** — field states cannot escape the ternary alphabet under the dynamics.

The leapfrog (Störmer-Verlet) integrator used here is:

```
E(t + dt/2) = E(t) + (dt/2) * curl B(t)
B(t + dt)   = B(t) - dt     * curl E(t + dt/2)
E(t + dt)   = E(t + dt/2) + (dt/2) * curl B(t + dt)
```

This is time-reversible and conserves the discrete EM energy to machine precision.

---

## API Overview

### `ElectricField`

Ternary charge distribution with Coulomb interactions. Each charge must be in {-1, 0, +1}.

```rust
use ternary_electromagnetism::ElectricField;

// Create a dipole: +1 charge at (0,0), -1 charge at (2,0)
let ef = ElectricField::new(
    vec![1, -1],
    vec![(0.0, 0.0), (2.0, 0.0)],
);

// Coulomb force between two unit charges at distance 1.0
// F = k*q1*q2/r^2 = 1.0*1*1/1.0 = 1.0 (repulsive)
let f = ElectricField::coulomb_force(1, 1, 1.0);
assert!(f > 0.0);

// Electric field vector at point (1.0, 1.0)
let (ex, ey) = ef.field_at((1.0, 1.0));
println!("E = ({:.4}, {:.4})", ex, ey);
```

### `MagneticField`

Ternary current sources using the Biot-Savart law. Each current value must be in {-1, 0, +1}.

```rust
use ternary_electromagnetism::MagneticField;

// Single wire at origin, current = +1 (out of page)
let mf = MagneticField::new(vec![1], vec![(0.0, 0.0)]);

// B at distance r=1: B = mu*I/(2*pi*r) = 1/(2*pi) ≈ 0.1592
let bz = mf.field_at((1.0, 0.0));
println!("Bz = {:.6}", bz); // 0.159155...

// Two wires with opposite currents: fields cancel at midpoint
let mf2 = MagneticField::new(vec![1, -1], vec![(0.0, 0.0), (2.0, 0.0)]);
let midpoint_bz = mf2.field_at((1.0, 0.0));
assert!(midpoint_bz.abs() < 1e-10); // ≈ 0
```

### `YeeLattice`

Discrete Maxwell's equations on an N×N Yee lattice with Ex, Ey, and Bz components.

```rust
use ternary_electromagnetism::YeeLattice;

// Create a 16x16 lattice
let mut lat = YeeLattice::new(16);

// Inject a field perturbation
lat.ex[8][8] = 1.0;

// Individual updates (for custom stepping)
lat.update_e(0.05);  // update E from curl B
lat.update_b(0.1);   // update B from curl E
lat.update_e(0.05);  // half-step to complete leapfrog

// Or use the built-in leapfrog step
let mut lat2 = YeeLattice::new(16);
lat2.ex[8][8] = 1.0;
lat2.step(0.1);  // full leapfrog: E/2 → B → E/2
```

### `WavePropagation`

High-level EM wave propagation through a ternary Yee lattice.

```rust
use ternary_electromagnetism::WavePropagation;

// Create propagation on a 32x32 lattice
let mut wp = WavePropagation::new(32);

// Check initial energy (should be ~0)
println!("Initial energy: {:.2e}", wp.energy()); // ~0

// Inject a pulse at cell (16, 16)
wp.inject_pulse(16, 16, 1.0);
println!("Energy after pulse: {:.4}", wp.energy()); // 0.5

// Advance 100 steps with dt=0.1
wp.advance(100, 0.1);
println!("Steps completed: {}", wp.steps);     // 100
println!("Energy at t=10:  {:.4}", wp.energy()); // conserved (dispersed across lattice)
```

### `Polarization`

Ternary polarization states mapped to {-1, 0, +1} with Jones vector representation.

```rust
use ternary_electromagnetism::Polarization;

// Construct from ternary value
let h = Polarization::from_value(-1); // Horizontal
let n = Polarization::from_value(0);  // None
let v = Polarization::from_value(1);  // Vertical

// Round-trip
assert_eq!(h.to_value(), -1);
assert_eq!(v.to_value(),  1);

// Jones vectors
assert_eq!(h.jones_vector(), (1.0, 0.0)); // Ex = 1, Ey = 0
assert_eq!(v.jones_vector(), (0.0, 1.0)); // Ex = 0, Ey = 1
assert_eq!(n.jones_vector(), (0.0, 0.0)); // no field

// Pattern match
match Polarization::from_value(1) {
    Polarization::Vertical   => println!("vertically polarized"),
    Polarization::Horizontal => println!("horizontally polarized"),
    Polarization::None       => println!("unpolarized"),
}
```

### `Interference`

Double-slit experiment and ternary phase interference.

```rust
use ternary_electromagnetism::Interference;

// Double-slit: maximum at center (x=0)
let i_center = Interference::double_slit_intensity(0.0, 1.0, 1.0);
assert!((i_center - 1.0).abs() < 1e-10); // = 1.0

// First minimum at x = λ/(2d) = 0.5
let i_min = Interference::double_slit_intensity(0.5, 1.0, 1.0);
assert!(i_min.abs() < 1e-10); // = 0.0

// Ternary phase interference
let constructive = Interference::ternary_phase_interference(1, 0);  // +1
let destructive  = Interference::ternary_phase_interference(1, -1); // 0
let clamped      = Interference::ternary_phase_interference(1, 1);  // +2 → +1 (clamped)
println!("{} {} {}", constructive, destructive, clamped); // 1 0 1
```

---

## Double-Slit Example

The double-slit formula implemented here is the standard scalar intensity:

```
I(x) = cos²(π · d · x / λ)
```

where `d` is the slit separation and `λ` is the wavelength. The table below shows the pattern for `d = 1.0`, `λ = 1.0`:

| x    | π·d·x/λ | I(x)   | Description          |
|------|---------|--------|----------------------|
| 0.00 | 0       | 1.000  | Central maximum      |
| 0.25 | π/4     | 0.500  | Half-power point     |
| 0.50 | π/2     | 0.000  | First minimum        |
| 0.75 | 3π/4    | 0.500  | Half-power point     |
| 1.00 | π       | 1.000  | Second maximum       |
| 1.50 | 3π/2    | 0.000  | Second minimum       |

Running this in code:

```rust
use ternary_electromagnetism::Interference;

let d = 1.0_f64;
let lambda = 1.0_f64;

for i in 0..=8 {
    let x = i as f64 * 0.25;
    let intensity = Interference::double_slit_intensity(x, d, lambda);
    println!("x = {:.2}  I = {:.4}", x, intensity);
}
```

Expected output:
```
x = 0.00  I = 1.0000
x = 0.25  I = 0.5000
x = 0.50  I = 0.0000
x = 0.75  I = 0.5000
x = 1.00  I = 1.0000
x = 1.25  I = 0.5000
x = 1.50  I = 0.0000
x = 1.75  I = 0.5000
x = 2.00  I = 1.0000
```

---

## Running the Tests

```bash
# Run all 30 tests
cargo test

# Run a specific test
cargo test double_slit_intensity_maximum_at_center

# Run with output shown
cargo test -- --nocapture
```

The test suite covers:

| Category        | Tests |
|-----------------|-------|
| ElectricField   | charge validation, Coulomb sign/magnitude/singularity, field symmetry |
| MagneticField   | Biot-Savart singularity/sign/value, field summation cancellation |
| YeeLattice      | zero init, update_e/update_b correctness, full step |
| WavePropagation | zero initial energy, pulse injection, step counting, vacuum conservation |
| Polarization    | round-trip, Jones vectors for all 3 states, invalid value panic |
| Interference    | maximum at x=0, first minimum, range bounds, additive/clamp cases |

---

## Design Notes

### Why Ternary?

Ternary (balanced base-3) arithmetic has been studied since Setun (1958). In the EM context, restricting field values to {-1, 0, +1} provides:

1. **Exact arithmetic** — no floating-point drift in the field states themselves
2. **Symmetry** — the positive/negative symmetry of Maxwell's equations is preserved by construction
3. **Minimal encoding** — a single trit encodes field direction and presence simultaneously
4. **Hardware efficiency** — ternary logic circuits have theoretical advantages in wire count and power

The continuous magnitudes used internally (`f64` for positions, distances, and intermediate values) represent the underlying physical coordinates, not the field values themselves.

### Boundary Conditions

The `YeeLattice` uses zero-padding boundary conditions: update loops skip the outermost row/column, leaving those cells at zero. This is equivalent to placing perfectly absorbing boundaries at the lattice edges. A future version could support periodic (toroidal) or Mur absorbing boundary conditions.

### Stability Criterion

For the Yee scheme in 2D, the Courant-Friedrichs-Lewy (CFL) stability criterion requires:

```
dt ≤ dx / (c * sqrt(2))
```

For the default lattice with unit cell spacing (dx = 1) and c = 1, use `dt ≤ 0.707`. The examples use `dt = 0.1` for safety.

---

## License

MIT
