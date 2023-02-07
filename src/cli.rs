use clap::Parser;

#[derive(Debug, Parser)]
#[command(version, author, about)]
pub struct Cli {
    #[arg(short = 'c', long)]
    pub configuration: String,
}
