use crate::consts::*;
use crate::vec3::Vec3;

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
        for i in 0..self.positions.len() {
            for j in (i + 1)..self.positions.len() {
                let p_i = self.positions[i];
                let p_j = self.positions[j];

                let dist = R_STAR.powi(2) / p_i.dist_sq(p_j);
                self.energy += dist.powi(6) - 2.0 * dist.powi(3);

                let du_ij = Vec3::broadcast(-48.0 * EPSILON_STAR * (dist.powi(7) - dist.powi(4)));
                let force = du_ij * (p_i - p_j);
                self.forces[i] += force;
                self.forces[j] += -force;
            }
        }

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
            for i in 0..self.positions.len() {
                for j in 0..self.positions.len() {
                    if i == j {
                        continue;
                    }

                    let tmp = self.positions[j] + trans_vec[k];
                    let dist = self.positions[i].dist_sq(tmp);
                    if dist > R_CUT.powi(2) {
                        continue;
                    }
                    let dist = R_STAR.powi(2) / dist;
                    self.energy += dist.powi(6) - 2.0 * dist.powi(3);

                    let du_ij =
                        Vec3::broadcast(-48.0 * EPSILON_STAR * (dist.powi(7) - dist.powi(4)));

                    self.forces[i] += du_ij * (self.positions[i] - tmp);
                }
            }
        }

        self.sum_of_forces = self.forces.iter().fold(Vec3::zero(), |acc, f| acc + *f);
        self.energy *= 2.0 * EPSILON_STAR;
    }
}
