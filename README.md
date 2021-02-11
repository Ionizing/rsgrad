![Rust](https://github.com/Ionizing/rsgrad/workflows/Rust/badge.svg)

Tracking the relaxation or MD progress of VASP calculation.

# Usage

```
$ ./rsgrad
rsgrad 0.2
Ionizing PeterSmith_9@outlook.com
Tracking the relaxation or MD progress of VASP calculation

USAGE:
    rsgrad [OPTIONS] [input]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --print-energy <t|f>        Print TOTEN in eV [default: f]
        --print-energyz <t|f>       Print TOTEN without entropy in eV [default: t]
        --print-favg <t|f>          Print averaged total force in eV/A [default: t]
        --print-fmax <t|f>          Print maximum total force in eV/A [default: t]
        --print-fmax-axis <t|f>     Print the axis where the strongest force weight lies on [default: f]
        --print-fmax-index <t|f>    Print the index of ion with maximum force load [default: f]
        --print-lgde <t|f>          Print Log10(deltaE), where deltaE is the
                                        difference between two consecutive energies
                                        without entropy [default: t]
        --print-magmom <t|f>        Print the magnetic moment in current system [default: t]
        --print-nscf <t|f>          Print the number SCF iteration in the ionic iteration [default: t]
        --print-time <t|f>          Print the time usage of each ionic iteration in unit of minute [default: t]
        --print-volume <t|f>        Print the system volumes of each ionic iteration [default: f]

ARGS:
    <input>     [default: OUTCAR]
```

# Origin

This tool is originated from the previous [repository](https://github.com/Ionizing/usefultools-for-vasp).

Now it is rewritten with Rust language, which make it safer and faster and more importantly, easier to develop.

# Ability

- Display the TOTEN and TOTEN without entropy info
- Display the number of SCF steps
- Display magnetic moment info
- Display TOTAL-FORCE info, including averag force and maximum force
- Display the time usage of each ionic step

# Future features
- [ ] A prettier output layout
- [ ] Output the trajectory of relaxation or MD
- [ ] Display the unconverged atoms
- [ ] Save the viberation modes
- [ ] More detailed error messages

# How to build

- Make sure the public network is accesible to your machine
- If you still haven't installed Rust on your machine, follow [this](https://www.rust-lang.org/tools/install) to install
- Clone this repo, and `cd rsgrad`
- To run all tests: `cargo test`
- To get the binay: `cargo build --release`, and the `rsgrad` is on `rsgrad/target/release/rsgrad`
- Have fun!
