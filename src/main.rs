pub mod cli;
pub mod consts;
pub mod particle_system;
pub mod vec3;

use crate::cli::{AlgorithmKind, Cli};
use crate::particle_system::ParticleSystem;

use clap::Parser;
use std::time::Instant;

fn main() {
    let args = Cli::parse();

    let mut system = ParticleSystem::from_config(
        args.position_configuration.into(),
        args.kinetic_momentums_configuration.into(),
    );

    match args.algorithm {
        AlgorithmKind::ClassicalLennardJones => {
            let now = Instant::now();
            system.lennard_jones();
            let duration = now.elapsed();
            println!(
                "Lennard-Jones (finished in {:?}):\n\tPotential energy = {:.5} eV\n\t∑ forces: {:?} N",
                duration, system.potential_energy, system.sum_of_forces
            );
        }
        AlgorithmKind::PeriodicalLennardJones => {
            let now = Instant::now();
            system.periodical_lennard_jones();
            let duration = now.elapsed();
            println!(
                "Periodical Lennard-Jones (finished in {:?}):\n\tPotential energy = {:.5} eV\n\t∑ forces: {:?} N",
                duration, system.potential_energy, system.sum_of_forces
            );
        }
        AlgorithmKind::ClassicalVerletVelocity => {
            // Initializes kinetic momentum for all particles of the system.
            system.init_kinetic_momentums();
            // Compute forces using periodical Lennard-Jones algorithm.
            system.lennard_jones();
            system.kinetic_energy = system.compute_kinetic_energy();
            system.temperature = system.compute_temperature();

            let now = Instant::now();
            for i in 0..args.n_step {
                system.verlet_velocity();

                system.lennard_jones();
                system.kinetic_energy = system.compute_kinetic_energy();
                system.temperature = system.compute_temperature();

                if i % args.m_step == 0 {
                    system.berendsen_thermostat();
                }

                println!(
                    "  Step {}, E = {:.5} eV, ∑f: {:.5e} N, KE: {:.5} J, T: {:.2} K",
                    i,
                    system.potential_energy,
                    system.sum_of_forces.sum(),
                    system.kinetic_energy,
                    system.temperature
                );

                system.save_simulation_iteration(args.output.as_str(), i);
            }
            let duration = now.elapsed();
            println!(
                "Verlet velocity (finished in {:?}):\n\tE= {:.5} eV\n\t∑ forces: {:.5e} N\n\tKinetic energy: {} J\n\tTemperature: {:.2} K",
                duration, system.potential_energy, system.sum_of_forces.sum(), system.kinetic_energy, system.temperature
            );
        }
        AlgorithmKind::PeriodicalVerletVelocity => {
            // system.init_kinetic_momentums();
            system.kinetic_energy = system.compute_kinetic_energy();
            system.temperature = system.compute_temperature();

            // Compute forces using periodical Lennard-Jones algorithm.
            system.periodical_lennard_jones();

            // CSV header
            println!("Step; Total energy; Potential energy; Temperature;");
            let now = Instant::now();
            for i in 0..args.n_step {
                system.verlet_velocity();
                system.kinetic_energy = system.compute_kinetic_energy();
                system.temperature = system.compute_temperature();

                system.periodical_lennard_jones();

                if i % args.m_step == 0 {
                    system.berendsen_thermostat();
                }

                println!(
                    "{}; {}; {}; {};",
                    i,
                    system.potential_energy + system.kinetic_energy,
                    system.potential_energy,
                    system.temperature
                );

                system.save_simulation_iteration(args.output.as_str(), i);
            }
            let duration = now.elapsed();
            println!(
                "Verlet velocity (finished in {:?}):\n\tE= {:.5} eV\n\t∑ forces: {:.5e} N\n\tKinetic energy: {} J\n\tTemperature: {:.2} K",
                duration, system.potential_energy, system.sum_of_forces.sum(), system.kinetic_energy, system.temperature
            );
        }
        AlgorithmKind::Exam => {
            // println!("rc; m;");

            // let rc = 8.0;
            // let m = system.mean_nb_particles_in_range(rc);
            // println!("{rc}; {m};");

            // let rc = 10.0;
            // let m = system.mean_nb_particles_in_range(rc);
            // println!("{rc}; {m};");

            // let rc = 12.0;
            // let m = system.mean_nb_particles_in_range(rc);
            // println!("{rc}; {m};");

            // let rc = 14.0;
            // let m = system.mean_nb_particles_in_range(rc);
            // println!("{rc}; {m};");

            // let sum_of_velocities = system.sum_of_velocities();
            // println!(
            //     "\n∑ velocities:\n\tx = {:.6e}\n\ty = {:.6e}\n\tz = {:.6e}",
            //     sum_of_velocities.x, sum_of_velocities.y, sum_of_velocities.z
            // );

            // // update kinetic energy
            // system.kinetic_energy = system.compute_kinetic_energy();
            // system.temperature = system.compute_temperature();
            // println!("temperature: t = {} k", system.temperature);

            // Compute forces using periodical Lennard-Jones algorithm.
            system.lennard_jones();

            // CSV header
            println!("Step; Total energy; Potential energy; Temperature;");
            let now = Instant::now();
            for i in 0..args.n_step {
                system.verlet_velocity();
                system.kinetic_energy = system.compute_kinetic_energy();
                system.temperature = system.compute_temperature();
                system.lennard_jones();

                println!(
                    "{}; {}; {}; {};",
                    i,
                    system.potential_energy + system.kinetic_energy,
                    system.potential_energy,
                    system.temperature
                );
            }
            let duration = now.elapsed();
            eprintln!("Simulation finished in {duration:?}");
        }
    }
}
