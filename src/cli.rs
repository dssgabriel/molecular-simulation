use std::{fmt::Display, str::FromStr};

use clap::Parser;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AlgorithmKind {
    ClassicalLennardJones,
    PeriodicalLennardJones,
    ClassicalVerletVelocity,
    PeriodicalVerletVelocity,
    Exam,
}

impl FromStr for AlgorithmKind {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "classical-lennard-jones" | "clj" => Ok(AlgorithmKind::ClassicalLennardJones),
            "periodical-lennard-jones" | "plj" => Ok(AlgorithmKind::PeriodicalLennardJones),
            "classical-verlet-velocity" | "cvv" => Ok(AlgorithmKind::ClassicalVerletVelocity),
            "periodical-verlet-velocity" | "pvv" => Ok(AlgorithmKind::PeriodicalVerletVelocity),
            "exam" | "e" => Ok(AlgorithmKind::Exam),
            _ => Err("Unknown algorithm kind"),
        }
    }
}

impl Display for AlgorithmKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ClassicalLennardJones => write!(f, "classical-lennard-jones"),
            Self::PeriodicalLennardJones => write!(f, "periodical-lennard-jones"),
            Self::ClassicalVerletVelocity => {
                write!(f, "classical-verlet-velocity")
            }
            Self::PeriodicalVerletVelocity => {
                write!(f, "periodical-verlet-velocity")
            }
            Self::Exam => write!(f, "exam"),
        }
    }
}

#[derive(Debug, Parser)]
#[command(version, author, about)]
pub struct Cli {
    /// Path to configuration file with particles positions.
    #[arg(short, long, default_value_t = String::from("exam270323.xyz"))]
    pub position_configuration: String,

    /// Path to configuration file with particles kinetic momentums.
    #[arg(short, long, default_value_t = String::from("exam270323.mci"))]
    pub kinetic_momentums_configuration: String,

    /// Path to output PDB file.
    #[arg(short, long, default_value_t = String::from("simulation.pdb"))]
    pub output: String,

    /// Number of simulation steps.
    #[arg(short, long, default_value_t = 10_000)]
    pub n_step: usize,

    /// Number of simulation steps in-between Berendsen thermostat.
    #[arg(short, long, default_value_t = 100)]
    pub m_step: usize,

    /// Algorithm to run.
    #[arg(short, long, default_value_t = AlgorithmKind::Exam)]
    pub algorithm: AlgorithmKind,
}
