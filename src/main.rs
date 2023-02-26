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
    // let mut system2 = system.clone();
    let mut system3 = system.clone();

    // Lennard-Jones
    let now = Instant::now();
    system.lennard_jones();
    let duration = now.elapsed();
    println!(
        "Lennard-Jones (finished in {:?}):\n\tEnergy = {:.5} eV\n\tSum of forces: {}",
        duration, system.energy, system.sum_of_forces
    );

    // Periodical Lennard-Jones
    // let now = Instant::now();
    // system2.periodical_lennard_jones();
    // let duration = now.elapsed();
    // println!(
    //     "\nPeriodical Lennard-Jones (finished in {:?}):\n\tEnergy = {:.5} eV\n\tSum of forces: {}",
    //     duration, system2.energy, system2.sum_of_forces
    // );

    // Verlet velocity
    // Initializes kinetic momentum for all particles of the system.
    system3.init_kinetic_momentum();

    // Compute forces using periodical Lennard-Jones algorithm.
    system3.periodical_lennard_jones();

    system3.update_kinetic_energy();
    system3.update_temperature();

    let now = Instant::now();
    for i in 0..10_000 {
        system3.verlet_velocity();
        system3.update_kinetic_energy();
        system3.update_temperature();
        // println!(
        //     "\tStep {}, energy = {:.5} eV, sum of forces: {}, kinetic energy: {}, temperature: {}",
        //     i, system3.energy, system3.sum_of_forces, system3.kinetic_energy, system3.temperature
        // );

        if i % 100 == 0 {
            system3.berendsen_thermostat();
        }
    }
    let duration = now.elapsed();
    println!(
        "\nVelocity verlet (finished in {:?}):\n\tEnergy = {:.5} eV\n\tSum of forces: {}\n\tKinetic energy: {}\n\tTemperature: {}",
        duration, system3.energy, system3.sum_of_forces, system3.kinetic_energy, system3.temperature
    );
}
