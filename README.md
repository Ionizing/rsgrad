Tracking the relaxation or MD progress of VASP calculation.

# Usage

```
$ ./rsgrad
USAGE:
    rsgrad [input]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

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
