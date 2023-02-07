pub mod cli;
pub mod consts;
pub mod particle_system;
pub mod vec3;

use crate::cli::Cli;
use crate::particle_system::ParticleSystem;

use clap::Parser;
use std::time::Instant;

fn main() {
    let args = Cli::parse();
    let mut system = ParticleSystem::from_config(args.configuration.into());
    let mut system2 = system.clone();

    // Lennard-Jones
    let now = Instant::now();
    system.lennard_jones();
    let duration = now.elapsed();
    println!(
        "Lennard-Jones (finished in {:?}):\n\tEnergy = {:.5} eV\n\tSum of forces: {}",
        duration, system.energy, system.sum_of_forces
    );

    // TODO: Periodical Lennard-Jones
    let now = Instant::now();
    system2.periodical_lennard_jones();
    let duration = now.elapsed();
    println!(
        "\nPeriodical Lennard-Jones (finished in {:?}):\n\tEnergy = {:.5} eV\n\tSum of forces: {}",
        duration, system2.energy, system2.sum_of_forces
    );
}
