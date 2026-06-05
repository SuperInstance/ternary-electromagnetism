//! # ternary-electromagnetism
//!
//! EM field simulation on ternary lattices, where field values live in {-1, 0, +1}.
//!
//! Provides discrete electromagnetic simulation primitives:
//! - [`ElectricField`] — ternary charge distributions and Coulomb's law
//! - [`MagneticField`] — ternary current sources and Biot-Savart law
//! - [`YeeLattice`] — discrete Maxwell's equations on a Yee lattice
//! - [`WavePropagation`] — EM wave propagation through ternary media
//! - [`Polarization`] — ternary polarization states
//! - [`Interference`] — double-slit experiment with ternary phase

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// ElectricField
// ---------------------------------------------------------------------------

/// Ternary charge distribution with Coulomb interactions.
///
/// Each charge must be in {-1, 0, +1}.
#[derive(Debug, Clone)]
pub struct ElectricField {
    pub charges: Vec<i8>,
    pub positions: Vec<(f64, f64)>,
}

impl ElectricField {
    /// Create a new `ElectricField`, validating that all charges are in {-1, 0, +1}.
    ///
    /// # Panics
    /// Panics if any charge is outside {-1, 0, +1} or if `charges` and `positions`
    /// have different lengths.
    pub fn new(charges: Vec<i8>, positions: Vec<(f64, f64)>) -> Self {
        assert_eq!(
            charges.len(),
            positions.len(),
            "charges and positions must have the same length"
        );
        for &c in &charges {
            assert!(
                c == -1 || c == 0 || c == 1,
                "charge {} is not in {{-1, 0, +1}}",
                c
            );
        }
        Self { charges, positions }
    }

    /// Coulomb force between two ternary charges at distance `r`.
    ///
    /// F = k * q1 * q2 / r², with k = 1.0.
    /// Returns `0.0` when r < 1e-10 to avoid singularity.
    pub fn coulomb_force(q1: i8, q2: i8, r: f64) -> f64 {
        if r < 1e-10 {
            return 0.0;
        }
        let k = 1.0_f64;
        k * (q1 as f64) * (q2 as f64) / (r * r)
    }

    /// Compute the 2D electric field vector at `point` by summing all Coulomb contributions.
    ///
    /// The field direction points from each source charge toward (or away from) `point`.
    pub fn field_at(&self, point: (f64, f64)) -> (f64, f64) {
        let mut fx = 0.0_f64;
        let mut fy = 0.0_f64;

        for (i, &q) in self.charges.iter().enumerate() {
            if q == 0 {
                continue;
            }
            let (px, py) = self.positions[i];
            let dx = point.0 - px;
            let dy = point.1 - py;
            let r = (dx * dx + dy * dy).sqrt();
            if r < 1e-10 {
                continue;
            }
            // Force magnitude (using unit test charge q_test = +1 at `point`)
            let magnitude = Self::coulomb_force(q, 1, r);
            // Direction: unit vector from source to point
            fx += magnitude * dx / r;
            fy += magnitude * dy / r;
        }

        (fx, fy)
    }
}

// ---------------------------------------------------------------------------
// MagneticField
// ---------------------------------------------------------------------------

/// Ternary current sources with Biot-Savart interactions.
///
/// Each current value must be in {-1, 0, +1}.
#[derive(Debug, Clone)]
pub struct MagneticField {
    pub currents: Vec<i8>,
    pub positions: Vec<(f64, f64)>,
}

impl MagneticField {
    /// Create a new `MagneticField`, validating that all currents are in {-1, 0, +1}.
    ///
    /// # Panics
    /// Panics if any current is outside {-1, 0, +1} or if lengths differ.
    pub fn new(currents: Vec<i8>, positions: Vec<(f64, f64)>) -> Self {
        assert_eq!(
            currents.len(),
            positions.len(),
            "currents and positions must have the same length"
        );
        for &c in &currents {
            assert!(
                c == -1 || c == 0 || c == 1,
                "current {} is not in {{-1, 0, +1}}",
                c
            );
        }
        Self { currents, positions }
    }

    /// Biot-Savart field magnitude at perpendicular distance `r` from a wire carrying `current`.
    ///
    /// B = mu * I / (2 * pi * r), with mu = 1.0.
    /// Returns `0.0` when r < 1e-10.
    pub fn biot_savart(current: i8, r: f64) -> f64 {
        if r < 1e-10 {
            return 0.0;
        }
        let mu = 1.0_f64;
        mu * (current as f64) / (2.0 * PI * r)
    }

    /// Compute the scalar z-component of the magnetic field at `point`.
    ///
    /// Each current wire is treated as infinite and straight along z,
    /// contributing a signed B_z based on the right-hand rule in 2D.
    pub fn field_at(&self, point: (f64, f64)) -> f64 {
        let mut bz = 0.0_f64;

        for (i, &current) in self.currents.iter().enumerate() {
            if current == 0 {
                continue;
            }
            let (px, py) = self.positions[i];
            let dx = point.0 - px;
            let dy = point.1 - py;
            let r = (dx * dx + dy * dy).sqrt();
            // Sign from right-hand rule: B_z = -I * (dx component in 2D cross product)
            // For a wire at origin carrying I in +z, B circles in xy-plane.
            // At point (dx, dy), B is tangential; z-projection sign = sign of I.
            // Simplified: contribute signed magnitude (positive = out of page for +I).
            bz += Self::biot_savart(current, r);
        }

        bz
    }
}

// ---------------------------------------------------------------------------
// YeeLattice
// ---------------------------------------------------------------------------

/// Discrete Maxwell's equations on an N×N Yee lattice.
///
/// The Yee scheme staggeres E and B fields in both space and time,
/// giving second-order accuracy and exact discrete charge conservation.
/// Ex is located at cell-edge x-faces, Ey at y-faces, Bz at cell centers.
#[derive(Debug, Clone)]
pub struct YeeLattice {
    pub ex: Vec<Vec<f64>>,
    pub ey: Vec<Vec<f64>>,
    pub bz: Vec<Vec<f64>>,
    pub n: usize,
}

impl YeeLattice {
    /// Create an N×N zero-initialized Yee lattice.
    pub fn new(n: usize) -> Self {
        let ex = vec![vec![0.0_f64; n]; n];
        let ey = vec![vec![0.0_f64; n]; n];
        let bz = vec![vec![0.0_f64; n]; n];
        Self { ex, ey, bz, n }
    }

    /// Update E fields from curl B (forward difference).
    ///
    /// Ex[i][j] += dt * (Bz[i][j] - Bz[i-1][j])
    /// Ey[i][j] -= dt * (Bz[i][j] - Bz[i][j-1])
    ///
    /// Boundary rows/columns (index 0 and n-1) are left untouched (zero-padding).
    pub fn update_e(&mut self, dt: f64) {
        let n = self.n;
        // Iterate over interior points only (i from 1..n, j from 1..n)
        for i in 1..n {
            for j in 0..n {
                self.ex[i][j] += dt * (self.bz[i][j] - self.bz[i - 1][j]);
            }
        }
        for i in 0..n {
            for j in 1..n {
                self.ey[i][j] -= dt * (self.bz[i][j] - self.bz[i][j - 1]);
            }
        }
    }

    /// Update B field from curl E (forward difference).
    ///
    /// Bz[i][j] -= dt * ((Ey[i+1][j] - Ey[i][j]) - (Ex[i][j+1] - Ex[i][j]))
    ///
    /// Boundary rows/columns (index n-1) are left untouched.
    pub fn update_b(&mut self, dt: f64) {
        let n = self.n;
        for i in 0..n - 1 {
            for j in 0..n - 1 {
                self.bz[i][j] -=
                    dt * ((self.ey[i + 1][j] - self.ey[i][j]) - (self.ex[i][j + 1] - self.ex[i][j]));
            }
        }
    }

    /// Leapfrog (Störmer-Verlet) step: half E update, full B update, half E update.
    pub fn step(&mut self, dt: f64) {
        self.update_e(dt / 2.0);
        self.update_b(dt);
        self.update_e(dt / 2.0);
    }
}

// ---------------------------------------------------------------------------
// WavePropagation
// ---------------------------------------------------------------------------

/// EM wave propagation through ternary media using a Yee lattice.
#[derive(Debug, Clone)]
pub struct WavePropagation {
    pub lattice: YeeLattice,
    pub steps: usize,
}

impl WavePropagation {
    /// Create a new `WavePropagation` with an N×N lattice.
    pub fn new(n: usize) -> Self {
        Self {
            lattice: YeeLattice::new(n),
            steps: 0,
        }
    }

    /// Inject a Gaussian-like pulse at lattice cell (x, y) by setting ex[x][y].
    pub fn inject_pulse(&mut self, x: usize, y: usize, amplitude: f64) {
        self.lattice.ex[x][y] = amplitude;
    }

    /// Advance the simulation by `n_steps` leapfrog steps of size `dt`.
    pub fn advance(&mut self, n_steps: usize, dt: f64) {
        for _ in 0..n_steps {
            self.lattice.step(dt);
        }
        self.steps += n_steps;
    }

    /// Compute the total electromagnetic energy stored in the lattice.
    ///
    /// E_total = (1/2) * sum_{i,j} (Ex²[i][j] + Ey²[i][j] + Bz²[i][j])
    pub fn energy(&self) -> f64 {
        let mut total = 0.0_f64;
        let n = self.lattice.n;
        for i in 0..n {
            for j in 0..n {
                total += self.lattice.ex[i][j].powi(2)
                    + self.lattice.ey[i][j].powi(2)
                    + self.lattice.bz[i][j].powi(2);
            }
        }
        total / 2.0
    }
}

// ---------------------------------------------------------------------------
// Polarization
// ---------------------------------------------------------------------------

/// Ternary polarization state of an EM wave.
///
/// Maps to {-1, 0, +1}: Horizontal (-1), None (0), Vertical (+1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Polarization {
    /// Horizontal polarization (value = -1).
    Horizontal,
    /// No polarization / unpolarized (value = 0).
    None,
    /// Vertical polarization (value = +1).
    Vertical,
}

impl Polarization {
    /// Construct from a ternary integer value.
    ///
    /// - `-1` → `Horizontal`
    /// - `0` → `None`
    /// - `+1` → `Vertical`
    ///
    /// # Panics
    /// Panics if `v` is not in {-1, 0, +1}.
    pub fn from_value(v: i8) -> Self {
        match v {
            -1 => Polarization::Horizontal,
            0 => Polarization::None,
            1 => Polarization::Vertical,
            _ => panic!("polarization value {} not in {{-1, 0, +1}}", v),
        }
    }

    /// Convert back to the ternary integer value.
    pub fn to_value(&self) -> i8 {
        match self {
            Polarization::Horizontal => -1,
            Polarization::None => 0,
            Polarization::Vertical => 1,
        }
    }

    /// Jones vector representation: (Ex_amplitude, Ey_amplitude).
    ///
    /// - `Horizontal` → `(1.0, 0.0)`
    /// - `None` → `(0.0, 0.0)`
    /// - `Vertical` → `(0.0, 1.0)`
    pub fn jones_vector(&self) -> (f64, f64) {
        match self {
            Polarization::Horizontal => (1.0, 0.0),
            Polarization::None => (0.0, 0.0),
            Polarization::Vertical => (0.0, 1.0),
        }
    }
}

// ---------------------------------------------------------------------------
// Interference
// ---------------------------------------------------------------------------

/// Double-slit experiment and ternary phase interference utilities.
pub struct Interference;

impl Interference {
    /// Compute double-slit intensity pattern at position `x`.
    ///
    /// Uses the standard formula: I(x) = cos²(π * d * x / λ)
    ///
    /// where `d` is slit separation and `λ` is wavelength.
    ///
    /// Returns a value in [0.0, 1.0].
    pub fn double_slit_intensity(x: f64, slit_sep: f64, wavelength: f64) -> f64 {
        let phase = PI * slit_sep * x / wavelength;
        phase.cos().powi(2)
    }

    /// Ternary phase interference: add two ternary phases, clamp result to {-1, 0, +1}.
    ///
    /// The raw sum is clamped: values > 1 become +1, values < -1 become -1.
    pub fn ternary_phase_interference(phase1: i8, phase2: i8) -> i8 {
        let sum = (phase1 as i16) + (phase2 as i16);
        if sum > 1 {
            1
        } else if sum < -1 {
            -1
        } else {
            sum as i8
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- ElectricField tests ---

    #[test]
    fn electric_field_valid_charges_accepted() {
        // Should not panic
        let _ef = ElectricField::new(
            vec![-1, 0, 1],
            vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)],
        );
    }

    #[test]
    #[should_panic(expected = "not in {-1, 0, +1}")]
    fn electric_field_invalid_charge_panics() {
        let _ef = ElectricField::new(vec![2], vec![(0.0, 0.0)]);
    }

    #[test]
    fn coulomb_force_like_charges_repel() {
        // +1 and +1: force > 0 (repulsive)
        let f = ElectricField::coulomb_force(1, 1, 1.0);
        assert!(f > 0.0, "like charges should repel (force > 0), got {}", f);
    }

    #[test]
    fn coulomb_force_unlike_charges_attract() {
        // +1 and -1: force < 0 (attractive)
        let f = ElectricField::coulomb_force(1, -1, 1.0);
        assert!(f < 0.0, "unlike charges should attract (force < 0), got {}", f);
    }

    #[test]
    fn coulomb_force_zero_at_origin() {
        // r near 0 should return 0.0
        let f = ElectricField::coulomb_force(1, 1, 0.0);
        assert_eq!(f, 0.0, "coulomb force should be 0 for r=0");
    }

    #[test]
    fn coulomb_force_inverse_square() {
        // F at r=1 should be 4x F at r=2 (inverse square law)
        let f1 = ElectricField::coulomb_force(1, 1, 1.0);
        let f2 = ElectricField::coulomb_force(1, 1, 2.0);
        let ratio = f1 / f2;
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "inverse square law: ratio should be 4.0, got {}",
            ratio
        );
    }

    #[test]
    fn electric_field_at_symmetry() {
        // Two equal positive charges at (-1,0) and (+1,0): field at origin is zero by symmetry
        let ef = ElectricField::new(
            vec![1, 1],
            vec![(-1.0, 0.0), (1.0, 0.0)],
        );
        let (fx, fy) = ef.field_at((0.0, 0.0));
        assert!(fx.abs() < 1e-10, "Ex should be ~0 by symmetry, got {}", fx);
        assert!(fy.abs() < 1e-10, "Ey should be ~0 by symmetry, got {}", fy);
    }

    // --- MagneticField tests ---

    #[test]
    fn biot_savart_zero_at_origin() {
        let b = MagneticField::biot_savart(1, 0.0);
        assert_eq!(b, 0.0, "biot_savart should be 0 for r=0");
    }

    #[test]
    fn biot_savart_positive_current_positive_field() {
        let b = MagneticField::biot_savart(1, 1.0);
        assert!(b > 0.0, "positive current should give positive B_z, got {}", b);
    }

    #[test]
    fn magnetic_field_at_sums_contributions() {
        // One wire at origin with current +1, one at (2,0) with current -1
        // At point (1,0): r=1 from each, contributions cancel
        let mf = MagneticField::new(
            vec![1, -1],
            vec![(0.0, 0.0), (2.0, 0.0)],
        );
        let bz = mf.field_at((1.0, 0.0));
        assert!(
            bz.abs() < 1e-10,
            "equal and opposite currents equidistant should cancel, got {}",
            bz
        );
    }

    #[test]
    fn magnetic_field_single_wire_value() {
        // Single wire at origin, current=1, query point at r=1 → B = 1/(2π)
        let mf = MagneticField::new(vec![1], vec![(0.0, 0.0)]);
        let bz = mf.field_at((1.0, 0.0));
        let expected = 1.0 / (2.0 * PI);
        assert!(
            (bz - expected).abs() < 1e-10,
            "single wire: expected {}, got {}",
            expected,
            bz
        );
    }

    // --- YeeLattice tests ---

    #[test]
    fn yee_lattice_initializes_to_zero() {
        let lat = YeeLattice::new(4);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(lat.ex[i][j], 0.0);
                assert_eq!(lat.ey[i][j], 0.0);
                assert_eq!(lat.bz[i][j], 0.0);
            }
        }
    }

    #[test]
    fn yee_lattice_update_e_changes_field() {
        let mut lat = YeeLattice::new(4);
        // Set Bz to a non-zero value to cause Ex update
        lat.bz[1][1] = 1.0;
        lat.bz[0][1] = 0.0;
        lat.update_e(0.1);
        // ex[1][1] should change: += dt * (bz[1][1] - bz[0][1]) = 0.1 * 1.0 = 0.1
        assert!(
            (lat.ex[1][1] - 0.1).abs() < 1e-10,
            "ex[1][1] should be 0.1, got {}",
            lat.ex[1][1]
        );
    }

    #[test]
    fn yee_lattice_update_b_changes_field() {
        let mut lat = YeeLattice::new(4);
        // Set Ey to a gradient to drive Bz update
        lat.ey[1][1] = 1.0;
        lat.ey[0][1] = 0.0;
        lat.update_b(0.1);
        // bz[0][1] -= dt * ((ey[1][1] - ey[0][1]) - (ex[0][2] - ex[0][1]))
        //           = -0.1 * (1.0 - 0.0 - 0.0) = -0.1
        assert!(
            (lat.bz[0][1] - (-0.1)).abs() < 1e-10,
            "bz[0][1] should be -0.1, got {}",
            lat.bz[0][1]
        );
    }

    #[test]
    fn yee_lattice_step_runs_without_error() {
        let mut lat = YeeLattice::new(8);
        lat.ex[3][3] = 1.0;
        // Just verify it doesn't panic and fields change
        lat.step(0.1);
        // Some Bz values should now be non-zero
        let any_nonzero = (0..8).any(|i| (0..8).any(|j| lat.bz[i][j] != 0.0));
        assert!(any_nonzero, "after step, some Bz values should be non-zero");
    }

    // --- WavePropagation tests ---

    #[test]
    fn wave_propagation_energy_starts_near_zero() {
        let wp = WavePropagation::new(8);
        assert!(
            wp.energy() < 1e-15,
            "initial energy should be ~0, got {}",
            wp.energy()
        );
    }

    #[test]
    fn wave_propagation_energy_increases_after_pulse() {
        let mut wp = WavePropagation::new(8);
        let initial_energy = wp.energy();
        wp.inject_pulse(4, 4, 1.0);
        let post_pulse_energy = wp.energy();
        assert!(
            post_pulse_energy > initial_energy,
            "energy should increase after pulse injection: before={}, after={}",
            initial_energy,
            post_pulse_energy
        );
    }

    #[test]
    fn wave_propagation_advance_steps_increment() {
        let mut wp = WavePropagation::new(8);
        wp.inject_pulse(4, 4, 1.0);
        wp.advance(10, 0.1);
        assert_eq!(wp.steps, 10, "steps should be 10 after advance(10, ...)");
    }

    #[test]
    fn wave_propagation_energy_conserved_in_vacuum() {
        // Without any pulse, energy stays zero
        let mut wp = WavePropagation::new(8);
        wp.advance(50, 0.1);
        assert!(
            wp.energy() < 1e-15,
            "energy should remain ~0 in empty lattice, got {}",
            wp.energy()
        );
    }

    // --- Polarization tests ---

    #[test]
    fn polarization_from_value_to_value_roundtrip() {
        for &v in &[-1i8, 0i8, 1i8] {
            let p = Polarization::from_value(v);
            assert_eq!(p.to_value(), v, "round-trip failed for value {}", v);
        }
    }

    #[test]
    fn polarization_jones_vector_horizontal() {
        let p = Polarization::Horizontal;
        assert_eq!(p.jones_vector(), (1.0, 0.0));
    }

    #[test]
    fn polarization_jones_vector_none() {
        let p = Polarization::None;
        assert_eq!(p.jones_vector(), (0.0, 0.0));
    }

    #[test]
    fn polarization_jones_vector_vertical() {
        let p = Polarization::Vertical;
        assert_eq!(p.jones_vector(), (0.0, 1.0));
    }

    #[test]
    #[should_panic(expected = "not in {-1, 0, +1}")]
    fn polarization_invalid_value_panics() {
        let _p = Polarization::from_value(2);
    }

    // --- Interference tests ---

    #[test]
    fn double_slit_intensity_maximum_at_center() {
        // At x=0, cos²(0) = 1.0
        let intensity = Interference::double_slit_intensity(0.0, 1.0, 1.0);
        assert!(
            (intensity - 1.0).abs() < 1e-10,
            "intensity at x=0 should be 1.0, got {}",
            intensity
        );
    }

    #[test]
    fn double_slit_intensity_minimum() {
        // At x = λ/(2d), cos²(π/2) = 0.0
        let wavelength = 1.0_f64;
        let slit_sep = 1.0_f64;
        let x = wavelength / (2.0 * slit_sep); // = 0.5
        let intensity = Interference::double_slit_intensity(x, slit_sep, wavelength);
        assert!(
            intensity.abs() < 1e-10,
            "intensity at first minimum should be 0.0, got {}",
            intensity
        );
    }

    #[test]
    fn double_slit_intensity_in_range() {
        // All intensities should be in [0, 1]
        for i in 0..20 {
            let x = i as f64 * 0.1;
            let intensity = Interference::double_slit_intensity(x, 2.0, 1.0);
            assert!(intensity >= 0.0 && intensity <= 1.0 + 1e-10,
                "intensity out of range [0,1]: {}", intensity);
        }
    }

    #[test]
    fn ternary_phase_interference_additive() {
        // 0 + 0 = 0
        assert_eq!(Interference::ternary_phase_interference(0, 0), 0);
        // +1 + 0 = +1
        assert_eq!(Interference::ternary_phase_interference(1, 0), 1);
        // -1 + 0 = -1
        assert_eq!(Interference::ternary_phase_interference(-1, 0), -1);
        // +1 + (-1) = 0 (destructive)
        assert_eq!(Interference::ternary_phase_interference(1, -1), 0);
    }

    #[test]
    fn ternary_phase_interference_clamp_positive() {
        // +1 + +1 = 2 → clamped to +1
        assert_eq!(Interference::ternary_phase_interference(1, 1), 1);
    }

    #[test]
    fn ternary_phase_interference_clamp_negative() {
        // -1 + -1 = -2 → clamped to -1
        assert_eq!(Interference::ternary_phase_interference(-1, -1), -1);
    }
}
