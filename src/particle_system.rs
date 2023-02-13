use crate::consts::*;
use crate::vec3::Vec3;

use rayon::prelude::*;

use std::fs::File;
use std::io::Read;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParticleSystem {
    pub positions: Vec<Vec3>,
    pub forces: Vec<Vec3>,
    pub sum_of_forces: Vec3,
    pub energy: f64,
}

impl ParticleSystem {
    pub fn new() -> Self {
        ParticleSystem {
            positions: Vec::new(),
            forces: Vec::new(),
            sum_of_forces: Vec3::zero(),
            energy: 0.0,
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
}
