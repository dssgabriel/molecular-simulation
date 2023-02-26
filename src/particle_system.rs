use crate::consts::*;
use crate::vec3::Vec3;

use rand::distributions::{Distribution, Uniform};
use rayon::prelude::*;

use std::fs::File;
use std::io::Read;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParticleSystem {
    pub n_particles: usize,
    pub positions: Vec<Vec3>,
    pub forces: Vec<Vec3>,
    pub energy: f64,
    pub sum_of_forces: Vec3,
    pub kinetic_momentum: Vec<Vec3>,
    pub kinetic_energy: f64,
    pub temperature: f64,
}

impl ParticleSystem {
    pub fn new() -> Self {
        ParticleSystem {
            n_particles: 0,
            positions: Vec::new(),
            forces: Vec::new(),
            kinetic_momentum: Vec::new(),
            sum_of_forces: Vec3::zero(),
            energy: 0.0,
            kinetic_energy: 0.0,
            temperature: 0.0,
        }
    }

    pub fn from_config(filename: std::path::PathBuf) -> Self {
        let mut system = Self::new();
        let mut contents = String::new();

        let mut file = File::open(filename).expect("failed to open file `{filename}`");
        file.read_to_string(&mut contents)
            .expect("failed to read file `{filename}`");

        for line in contents.lines().skip(1) {
            let mut values = line
                .split_whitespace()
                .skip(1)
                .map(|val| val.parse::<f64>().expect("failed to parse {val}"));
            system.positions.push(Vec3::new(
                values.next().unwrap(),
                values.next().unwrap(),
                values.next().unwrap(),
            ));
            system.forces.push(Vec3::zero());
        }

        system
    }

    pub fn init_kinetic_momentum(&mut self) {
        let range = Uniform::from(-1.0_f64..1.0_f64);
        let mut rng = rand::thread_rng();

        // For each particle, generate a random kinetic momentum (P_x, P_y, P_z) in
        // the range of [-1.0, 1.0[.
        for _ in 0..self.n_particles {
            self.kinetic_momentum.push(Vec3::new(
                range.sample(&mut rng),
                range.sample(&mut rng),
                range.sample(&mut rng),
            ));
        }

        self.first_kinetic_momentum_recalibration();
        self.second_kinetic_momentum_recalibration();
        self.first_kinetic_momentum_recalibration();
    }

    fn first_kinetic_momentum_recalibration(&mut self) {
        // Recompute kinetic energy of the system
        self.update_kinetic_energy();

        let n_dl = 3.0 * self.n_particles as f64 - 3.0;
        let ratio = Vec3::splat((n_dl * R_CONST * T_0) / self.kinetic_energy);

        for p in &mut self.kinetic_momentum {
            *p *= ratio
        }
    }

    fn second_kinetic_momentum_recalibration(&mut self) {
        // Sum the kinetic momentum of the system on all axes.
        let mut sum_kinetic_momentum = Vec3::zero();
        for p in &self.kinetic_momentum {
            sum_kinetic_momentum += *p;
        }

        // Correct the kinetic momentum of the particles
        for p in &mut self.kinetic_momentum {
            *p -= sum_kinetic_momentum;
        }
    }

    pub fn update_kinetic_energy(&mut self) {
        for p in &self.positions {
            self.kinetic_energy += p.mag_sq() / M_I;
        }
        self.kinetic_energy *= 1.0 / 2.0 * CONVERSION_FORCE;
    }

    pub fn update_temperature(&mut self) {
        let n_dl = 3.0 * self.n_particles as f64 - 3.0;
        self.temperature = self.kinetic_energy * (1.0 / n_dl * R_CONST);
    }

    pub fn lennard_jones(&mut self) {
        self.energy += self
            .forces
            .par_iter_mut()
            .enumerate()
            .map(|(i, force_i)| {
                let mut energy = 0.0;

                for j in 0..self.positions.len() {
                    if i == j {
                        continue;
                    }

                    let distance =
                        R_STAR.powi(2) / self.positions[i].distance_square(self.positions[j]);
                    let du_ij =
                        Vec3::splat(-48.0 * EPSILON_STAR * (distance.powi(7) - distance.powi(4)));

                    // Compute force exerted on particle `i`
                    *force_i += du_ij * (self.positions[i] - self.positions[j]);

                    // Update energy of the particle system
                    // (but only for i < j so we don't double-count)
                    if i < j {
                        energy += distance.powi(6) - 2.0 * distance.powi(3);
                    }
                }

                energy
            })
            .sum::<f64>();

        self.sum_of_forces = self.forces.iter().fold(Vec3::zero(), |acc, f| acc + *f);
        self.energy *= 4.0 * EPSILON_STAR;
    }

    pub fn periodical_lennard_jones(&mut self) {
        // Initialize translation vector
        let mut trans_vec = Vec::with_capacity(N_SYM);
        for i in 0..trans_vec.capacity() {
            trans_vec.push(Vec3::new(
                (((i / 9) - 1) * 50) as f64,
                ((((i / 3) % 3) - 1) * 50) as f64,
                (((i % 3) - 1) * 50) as f64,
            ));
        }

        for k in 0..N_SYM {
            self.energy += self
                .forces
                .par_iter_mut()
                .enumerate()
                .map(|(i, force_i)| {
                    let mut energy = 0.0;

                    for j in 0..self.positions.len() {
                        if i == j {
                            continue;
                        }

                        let tmp = self.positions[j] + trans_vec[k];
                        let distance = self.positions[i].distance_square(tmp);
                        if distance > R_CUT.powi(2) {
                            continue;
                        }
                        let distance = R_STAR.powi(2) / distance;

                        let du_ij = Vec3::splat(
                            -48.0 * EPSILON_STAR * (distance.powi(7) - distance.powi(4)),
                        );

                        *force_i += du_ij * (self.positions[i] - tmp);
                        energy += distance.powi(6) - 2.0 * distance.powi(3);
                    }

                    energy
                })
                .sum::<f64>();
        }

        self.sum_of_forces = self.forces.iter().fold(Vec3::zero(), |acc, f| acc + *f);
        self.energy *= 2.0 * EPSILON_STAR;
    }

    pub fn verlet_velocity(&mut self) {
        // Initializes kinetic momentum for all particles of the system.
        self.init_kinetic_momentum();

        // Compute forces using periodical Lennard-Jones algorithm.
        self.periodical_lennard_jones();

        // Update the kinetic moments of each particle.
        let update = Vec3::splat(DT * CONVERSION_FORCE * 0.5);
        for (p, f) in self.kinetic_momentum.iter_mut().zip(self.forces.iter()) {
            *p -= update * *f;
        }

        // Update the positions of each particle.
        let vec_m_i = Vec3::splat(M_I);
        let vec_dt = Vec3::splat(DT);
        for (pos, p) in self.positions.iter_mut().zip(self.kinetic_momentum.iter()) {
            *pos += vec_dt * *p / vec_m_i;
        }

        // Re-compute forces.
        self.periodical_lennard_jones();

        // Update the kinetic moments of each particle, again.
        for (p, f) in self.kinetic_momentum.iter_mut().zip(self.forces.iter()) {
            *p -= update * *f;
        }
    }

    pub fn berendsen_thermostat(&mut self) {
        self.update_temperature();
        let vec_temp = Vec3::splat(self.temperature / T_0 - 1.0);
        let vec_gamma = Vec3::splat(GAMMA);

        for p in &mut self.kinetic_momentum {
            *p += *p * vec_gamma * vec_temp;
        }
    }
}
