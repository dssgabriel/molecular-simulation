#![allow(non_snake_case)]
use crate::consts::*;
use crate::vec3::Vec3;

use rand::distributions::{Distribution, Uniform};
use rayon::prelude::*;

use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::path::PathBuf;

/// A particle system for molecular dynamics.
///
/// The `positions`, `forces` and `kinetic_momentums` vectors are composed of `n_particles`
/// elements.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParticleSystem {
    pub n_particles: usize,
    pub positions: Vec<Vec3>,
    pub forces: Vec<Vec3>,
    pub kinetic_momentums: Vec<Vec3>,
    pub velocities: Vec<Vec3>,
    pub sum_of_forces: Vec3,
    pub potential_energy: f64,
    pub kinetic_energy: f64,
    pub temperature: f64,
}

impl ParticleSystem {
    /// Initializes the positions of particles in a particle system from a configuration file.
    pub fn from_config(filename_pos: PathBuf, filename_km: PathBuf) -> Self {
        let mut system = Self::default();
        let mut contents_pos = String::new();
        let mut contents_km = String::new();

        let mut file_pos = File::open(filename_pos).expect("failed to open file `{filename}`");
        file_pos
            .read_to_string(&mut contents_pos)
            .expect("failed to read file `{filename}`");

        let mut file_km = File::open(filename_km).expect("failed to open file `{filename}`");
        file_km
            .read_to_string(&mut contents_km)
            .expect("failed to read file `{filename}`");

        // Read each line of the and parse its values
        // Skip the first line (header) and the first element of each line (particle ID,
        // not used here)
        for (line_p, line_km) in contents_pos
            .lines()
            .skip(1)
            .zip(contents_km.lines().skip(1))
        {
            // Read positions of current line
            let mut pos = line_p
                .split_whitespace()
                .skip(1)
                .map(|val| val.parse::<f64>().expect("failed to parse {val}"));
            system.positions.push(Vec3::new(
                pos.next().unwrap(),
                pos.next().unwrap(),
                pos.next().unwrap(),
            ));

            // Read kinetic momentums of current line
            let mut km = line_km
                .split_whitespace()
                .skip(1)
                .map(|val| val.parse::<f64>().expect("failed to parse {val}"));
            system.kinetic_momentums.push(Vec3::new(
                km.next().unwrap(),
                km.next().unwrap(),
                km.next().unwrap(),
            ));

            // Compute velocity of particle on the current line given its kinetic momentum:
            //   km = v * m  =>  v = km / m
            let last_km = system.kinetic_momentums.last().unwrap();
            system.velocities.push(*last_km / M_I);

            // Initial force is set to zero
            system.forces.push(Vec3::zero());
        }

        system.n_particles = system.positions.len();

        // Make sure we have the same number of particles everywhere
        assert_eq!(system.n_particles, system.forces.len());
        assert_eq!(system.n_particles, system.kinetic_momentums.len());
        assert_eq!(system.n_particles, system.velocities.len());

        system
    }

    /// Initializes the kinetic momentum of each particle in the system.
    pub fn init_kinetic_momentums(&mut self) {
        // Kinetic momentums shall be in the range [-1.0, 1.0]
        let range = Uniform::from(-1.0_f64..1.0_f64);
        let mut rng = rand::thread_rng();

        // For each particle, generate a random kinetic momentum (P_x, P_y, P_z)
        for P in self.kinetic_momentums.iter_mut() {
            *P = Vec3::new(
                range.sample(&mut rng),
                range.sample(&mut rng),
                range.sample(&mut rng),
            );
        }

        // Perform the appropriate recalibrations
        self.first_kinetic_momentum_recalibration();
        self.second_kinetic_momentum_recalibration();
    }

    fn first_kinetic_momentum_recalibration(&mut self) {
        // Recompute kinetic energy of the system
        self.kinetic_energy = self.compute_kinetic_energy();

        let N_dl = (3 * self.n_particles - 3) as f64;
        let kinetic_energy_0 = N_dl * K_BOLTZMANN * T_0;
        let ratio = (kinetic_energy_0 / self.kinetic_energy).sqrt();

        for P in self.kinetic_momentums.iter_mut() {
            *P *= ratio
        }
    }

    fn second_kinetic_momentum_recalibration(&mut self) {
        // Sum the kinetic momentum of the system on all axes.
        let mut sum_kinetic_momentums = Vec3::zero();
        for P in self.kinetic_momentums.iter() {
            sum_kinetic_momentums += *P;
        }

        // Correct the kinetic momentum of the particles
        for P in self.kinetic_momentums.iter_mut() {
            *P -= sum_kinetic_momentums;
        }
    }

    /// Computes the kinetic energy of the system following the equation:
    ///   (1 / 2 * CONVERSION_FORCE) * ∑_i=1^N ((p_x^i)² + (p_y^i)² + (p_z^i)²) / m_i
    pub fn compute_kinetic_energy(&self) -> f64 {
        let mut sum = 0.0;
        for P in self.kinetic_momentums.iter() {
            sum += P.mag_sq() / M_I;
        }
        (CONVERSION_FORCE_2).recip() * sum
    }

    pub fn compute_temperature(&self) -> f64 {
        let N_dl = (3 * self.n_particles - 3) as f64;
        (N_dl * K_BOLTZMANN).recip() * self.kinetic_energy
    }

    pub fn mean_nb_particles_in_range(&self, cut_radius: f64) -> f64 {
        let mut count = 0;
        let cut_radius_sq = cut_radius.powi(2);

        for &pi in self.positions.iter() {
            for &pj in self.positions.iter() {
                if pi.distance_square(pj) < cut_radius_sq {
                    count += 1;
                }
            }
        }

        count as f64 / self.n_particles as f64
    }

    pub fn sum_of_velocities(&self) -> Vec3 {
        let mut sum = Vec3::zero();
        for &v in self.velocities.iter() {
            sum += v;
        }
        sum
    }

    pub fn lennard_jones(&mut self) {
        self.forces.par_iter_mut().for_each(|f| {
            *f = Vec3::zero();
        });
        self.sum_of_forces = Vec3::zero();

        self.potential_energy += self
            .forces
            .par_iter_mut()
            .enumerate()
            .map(|(i, force_i)| {
                let mut energy = 0.0;

                for j in 0..self.positions.len() {
                    if i == j {
                        continue;
                    }

                    let distance = R_STAR_SQ / self.positions[i].distance_square(self.positions[j]);
                    let distance_cube = distance.powi(3);
                    let du_ij = -48.0 * EPSILON_STAR * (distance.powi(4) * (distance_cube - 1.0));

                    // Compute force exerted on particle `i`
                    *force_i += (self.positions[i] - self.positions[j]) * du_ij;

                    // Update energy of the particle system
                    // (but only for i < j so we don't double-count)
                    if i < j {
                        energy += distance_cube * (distance_cube - 2.0);
                    }
                }

                energy
            })
            .sum::<f64>();

        self.sum_of_forces = self.forces.iter().fold(Vec3::zero(), |acc, f| acc + *f);
        self.potential_energy *= 4.0 * EPSILON_STAR;
    }

    pub fn periodical_lennard_jones(&mut self) {
        self.forces.par_iter_mut().for_each(|f| {
            *f = Vec3::zero();
        });
        self.sum_of_forces = Vec3::zero();

        for k in 0..N_SYM {
            self.potential_energy += self
                .forces
                .par_iter_mut()
                .enumerate()
                .map(|(i, force_i)| {
                    let mut energy = 0.0;

                    for j in 0..self.n_particles {
                        if i == j {
                            continue;
                        }

                        let translated_position_j = self.positions[j] + TRANSLATION_VECTORS[k];
                        let distance = self.positions[i].distance_square(translated_position_j);
                        if distance > R_CUT_SQ {
                            continue;
                        }
                        let distance = R_STAR_SQ / distance;

                        let distance_cube = distance.powi(3);
                        let du_ij =
                            -48.0 * EPSILON_STAR * (distance.powi(4) * (distance_cube - 1.0));

                        *force_i += (self.positions[i] - translated_position_j) * du_ij;

                        // Update energy of the particle system
                        // (but only for i < j so we don't double-count)
                        if i < j {
                            energy += distance_cube * (distance_cube - 2.0);
                        }
                    }

                    energy
                })
                .sum::<f64>();
        }

        self.sum_of_forces = self.forces.iter().fold(Vec3::zero(), |acc, f| acc + *f);
        self.potential_energy *= 4.0 * EPSILON_STAR;
    }

    pub fn verlet_velocity(&mut self) {
        let update = DT * CONVERSION_FORCE * 0.5;

        // Update the kinetic moments of each particle.
        for (P, f) in self.kinetic_momentums.iter_mut().zip(self.forces.iter()) {
            *P -= *f * update;
        }

        self.first_kinetic_momentum_recalibration();
        self.second_kinetic_momentum_recalibration();

        // Update the positions of each particle.
        for (pos, P) in self.positions.iter_mut().zip(self.kinetic_momentums.iter()) {
            *pos += (*P * DT) / M_I;
        }
    }

    pub fn berendsen_thermostat(&mut self) {
        for P in &mut self.kinetic_momentums {
            *P += GAMMA * (T_0 / self.temperature - 1.0);
        }
    }

    pub fn save_simulation_iteration(&self, filename: &str, iteration: usize) {
        let mut file = OpenOptions::new()
            .create(true)
            .write(true)
            .append(true)
            .open(filename)
            .unwrap();

        writeln!(file, "CRYST1  {L}  {L}  {L}  90.00  90.00  90.00  P  1").unwrap();
        writeln!(file, "MODEL  {iteration}").unwrap();
        for (i, p) in self.positions.iter().enumerate() {
            writeln!(file, "ATOM  {}  C  0  {} MRES", i + 1, *p).unwrap();
        }
        writeln!(file, "TER").unwrap();
        writeln!(file, "ENDMDL").unwrap();
    }
}
