//! Constants used for molecular dynamics simulations.

use crate::vec3::Vec3;

/// Lennard-Jones 𝜀* parameter.
pub const EPSILON_STAR: f64 = 0.2;

/// Lennard-Jones r* parameter.
pub const R_STAR: f64 = 3.0;
pub const R_STAR_SQ: f64 = R_STAR * R_STAR;

/// Cut radius to correct the potential energy of a particle.
pub const R_CUT: f64 = 10.0;
pub const R_CUT_SQ: f64 = R_CUT * R_CUT;

/// Dimensions of the replicated cubic box in periodcial conditions.
pub const L: f64 = 30.0;

/// Translation vectors for periodical conditions.
pub const TRANSLATION_VECTORS: [Vec3; 27] = [
    Vec3::zero(),
    Vec3::z(L as f64),
    Vec3::y(L as f64),
    Vec3::new(0.0, L as f64, L as f64),
    Vec3::x(L as f64),
    Vec3::new(L as f64, 0.0, L as f64),
    Vec3::new(L as f64, L as f64, 0.0),
    Vec3::splat(L as f64),
    Vec3::new(L as f64, L as f64, -L as f64),
    Vec3::new(L as f64, -L as f64, L as f64),
    Vec3::new(L as f64, -L as f64, -L as f64),
    Vec3::new(-L as f64, L as f64, L as f64),
    Vec3::new(-L as f64, L as f64, -L as f64),
    Vec3::new(-L as f64, -L as f64, L as f64),
    Vec3::splat(-L as f64),
    Vec3::z(-L as f64),
    Vec3::y(-L as f64),
    Vec3::new(0.0, -L as f64, -L as f64),
    Vec3::x(-L as f64),
    Vec3::new(-L as f64, 0.0, -L as f64),
    Vec3::new(-L as f64, -L as f64, 0.0),
    Vec3::new(L as f64, 0.0, -L as f64),
    Vec3::new(L as f64, -L as f64, 0.0),
    Vec3::new(-L as f64, 0.0, L as f64),
    Vec3::new(0.0, L as f64, -L as f64),
    Vec3::new(0.0, -L as f64, L as f64),
    Vec3::new(-L as f64, L as f64, 0.0),
];

/// Number of symmetries in periodical conditions.
pub const N_SYM: usize = 27;

/// Conversion factor to apply to units of force.
pub const CONVERSION_FORCE: f64 = 0.0004186;
pub const CONVERSION_FORCE_2: f64 = 2.0 * CONVERSION_FORCE;

/// Particle mass.
pub const M_I: f64 = 18.0;

/// Boltzman constant parameter used to compute the kinetic temperature of the particle system.
pub const K_BOLTZMANN: f64 = 8.31e-7;

/// Initial temperature of the system (in Kelvin).
pub const T_0: f64 = 300.0;

/// Constant parameter used in Berendsen thermostat.
pub const GAMMA: f64 = 0.01;

/// Time step in femtoseconds.
pub const DT: f64 = 1.0;
