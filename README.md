# Introduction to Molecular Simulation - Project

## Building
Make sure you have Rust and the `cargo` toolchain installed.
```sh
cargo build --release
```

## Executing
Execute with a given configuration file:
```sh
cargo run --release -- --configuration particules.xyz
# Or the shorthand:
cargo r -r -q -- -c particules.xyz
```
