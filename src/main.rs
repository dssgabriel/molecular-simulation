use clap::Parser;
use std::io::Read;

#[derive(Debug, Parser)]
#[command(version, author, about)]
struct Cli {
    #[arg(short = 'c', long)]
    configuration: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ParticleSystem {
    kind: Vec<i8>,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
}

impl ParticleSystem {
    pub fn new() -> Self {
        Self {
            kind: Vec::new(),
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
        }
    }

    pub fn from_config(filename: std::path::PathBuf) -> Self {
        let mut particle_sys = Self::new();
        let mut contents = String::new();

        let mut file = std::fs::File::open(filename).expect("failed to open file {filename}");
        file.read_to_string(&mut contents)
            .expect("failed to read file {filename}");

        for line in contents.lines().skip(1) {
            let mut line = line.split_whitespace();
            particle_sys.kind.push(
                line.next()
                    .unwrap()
                    .parse::<i8>()
                    .expect("failed to parse particle kind to `u8`"),
            );
            particle_sys.x.push(
                line.next()
                    .unwrap()
                    .parse::<f64>()
                    .expect("failed to parse particle's x coordinate to `f64`"),
            );
            particle_sys.y.push(
                line.next()
                    .unwrap()
                    .parse::<f64>()
                    .expect("failed to parse particle's y coordinate to `f64`"),
            );
            particle_sys.z.push(
                line.next()
                    .unwrap()
                    .parse::<f64>()
                    .expect("failed to parse particle's z coordinate to `f64`"),
            );
        }

        particle_sys
    }

    pub fn distance_squared(&self, idx1: usize, idx2: usize) -> f64 {
        let x1 = self.x[idx1];
        let x2 = self.x[idx2];
        let y1 = self.y[idx1];
        let y2 = self.y[idx2];
        let z1 = self.z[idx1];
        let z2 = self.z[idx2];

        (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1)
    }

    pub fn lennard_jones_potential(&self) -> f64 {
        let n = self.kind.len();
        let r_star_6 = 2.0 * 3.0_f64.powf(6.0);
        let r_star_12 = 3.0_f64.powf(12.0);
        let eps_star = 0.2_f64;

        let mut u_lj = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                // Compute r_ij distance squared
                let r_ij = 1.0 / self.distance_squared(i, j);

                let r_ij_6 = r_ij * r_ij * r_ij;
                let r_ij_12 = r_ij_6 * r_ij_6;

                u_lj += eps_star * (r_star_12 * r_ij_12) - (r_star_6 * r_ij_6);
            }
        }

        u_lj * 4.0
    }
}

fn main() {
    let args = Cli::parse();
    let sys = ParticleSystem::from_config(args.configuration.into());
    let u_lj = sys.lennard_jones_potential();
    println!("U_LJ = {:E} eV", u_lj);
}
