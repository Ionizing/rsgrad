![Rust](https://github.com/Ionizing/rsgrad/workflows/Rust/badge.svg)

Tracking the relaxation or MD progress of VASP calculation.

# NOTE: THIS PROJECT IS UNDER DEVELOPMENT, THE DOCUMENTS MAY BE OUT OF DATE.

# Download

Now the pre-built binaries are built only for `x86-64` architecture.

Head to https://github.com/Ionizing/rsgrad/releases/latest and download the pre-built binaries corresponding to your platform.

## For Linux users
- `linux-x86_64.tar.gz` is dynamically linked to `glibc` 2.18, with little higher performance. It cannot run on some old linux distributions (centos 7 or lower and so on)
- `linux-x86_64-musl.tar.gz` is full statically linked with `musl` libc. Though would be a little slower, it can be run on most linux distributions.

## For MacOS users
- `macos-x86_64.tar.gz` is for Intel chip. If your PC has M1 chip, you can run it with Rosetta as well. However, an ARM native binary is better, which means you need to build `rsgrad` from source.

The binary for M1 Mac is considered to be published in a few weeks.

# Usage （OUTDATED, NEEDS UPDATE）
## Main interface

```
rsgrad 0.2.5
@Ionizing github.com/Ionizing/rsgrad
A tool used to track the VASP calculation result

USAGE:
    rsgrad <SUBCOMMAND>

FLAGS:
    -h, --help
            Prints help information

    -V, --version
            Prints version information


SUBCOMMANDS:
    help    Prints this message or the help of the given subcommand(s)
    rlx     Tracking relaxation or MD progress
    traj    Operations about relaxation/MD trajectory
    vib     Tracking vibration information
```

## ~~List the brief info of current system~~

Available in the near future.

## See the relaxation info of OUTCAR

```
./rsgrad rlx -o OUTCAR_multiple_ionic_steps
[2022-02-06T12:19:58Z INFO  rsgrad::commands::rlx] Parsing file "OUTCAR_multiple_ionic_steps" and "./POSCAR"
  #Step  TOTEN_z/eV LgdE   Fmax #SCF Time/m Mag/muB
      1  -253.61859  2.4  0.000   14   2.45   NoMag
      2  -253.61023 -2.1  0.192    8   1.31   NoMag
      3  -253.61629 -2.2  0.501    7   1.34   NoMag
      4  -253.58960 -1.6  0.649    8   1.44   NoMag
      5  -253.64364 -1.3  0.001    7   1.29   NoMag
```

Full functions of `rsgrad rlx`

```
rsgrad-rlx 0.2.5
Tracking relaxation or MD progress.

Contains the evolution of energy, maximum of Hellmann-Feynman forces, magnetic moments and time usage of each ionic
step.

Hint: This command may require POSCAR for atom constraints information.

USAGE:
    rsgrad rlx [FLAGS] [OPTIONS]

FLAGS:
    -h, --help
            Prints help information

        --no-totenz
            Don't print TOTEN without entropy in eV
    
        --no-fmax
            Don't print maximum total force in A^3
    
        --no-lgde
            Don't print Log10(delta(TOTEN without entropy))
    
        --no-magmom
            Don't print total magnetic moment in muB
    
        --no-nscf
            Don't print number of SCF iteration for each ionic step
    
        --no-time
            Don't print time elapsed for each ionic step in minutes
    
    -e, --toten
            Prints TOTEN in eV
    
    -a, --favg
            Prints averaged total force in eV/A
    
    -x, --fmaxis
            Prints the axis where the strongest total force component lies on. [XYZ]
    
    -i, --fmidx
            Prints the index of ion with maximum total force load. Starts from 1
    
    -v, --volume
            Prints lattice volume in A^3
    
    -V, --version
            Prints version information


OPTIONS:
    -o, --outcar <outcar>
            Specify the input OUTCAR file [default: ./OUTCAR]

    -p, --poscar <poscar>
            Specify the input POSCAR file [default: ./POSCAR]

```


## Operations on vibration modes

```
$ ./rsgrad vib -o OUTCAR_vibrations -l
[2022-02-06T12:23:57Z INFO  rsgrad::commands::vib] Parsing file "OUTCAR_vibrations"
# --------------- Viberation modes for this system --------------- #
  ModeIndex:    1  Frequency/cm-1:    3627.91026  IsImagine: False
  ModeIndex:    2  Frequency/cm-1:    3620.67362  IsImagine: False
  ModeIndex:    3  Frequency/cm-1:    3431.76345  IsImagine: False
  ModeIndex:    4  Frequency/cm-1:    1551.74081  IsImagine: False
  ModeIndex:    5  Frequency/cm-1:    1537.18628  IsImagine: False
  ModeIndex:    6  Frequency/cm-1:     388.96334  IsImagine: False
  ModeIndex:    7  Frequency/cm-1:     370.87662  IsImagine: False
  ModeIndex:    8  Frequency/cm-1:     370.09082  IsImagine: False
  ModeIndex:    9  Frequency/cm-1:       0.65835  IsImagine: False
  ModeIndex:   10  Frequency/cm-1:       0.75226  IsImagine:  True
  ModeIndex:   11  Frequency/cm-1:       1.87333  IsImagine:  True
  ModeIndex:   12  Frequency/cm-1:     702.43818  IsImagine:  True
[2022-02-06T12:23:57Z INFO  rsgrad] Time used: 43.754564ms
```

Full usage

```
$ ./rsgrad vib --help
rsgrad-vib 0.2.5
Tracking vibration information.

For systems enabled vibration mode calculation, this command can extract phonon eigenvalues and phonon eigenvectors at
Gamma point.

USAGE:
    rsgrad vib [FLAGS] [OPTIONS]

FLAGS:
    -h, --help
            Prints help information

    -l, --list
            Shows vibration modes in brief
    
    -x, --save-as-xsfs
            Saves each selected modes to XSF file
    
    -V, --version
            Prints version information


OPTIONS:
    -o, --outcar <outcar>
            Specify the input OUTCAR file [default: ./OUTCAR]

        --save-in <save-in>
            Define where the files would be saved [default: .]
    
    -i, --select-indices <select-indices>...
            Selects the indices to operate.
    
            Step indices start from '1', if '0' is given, all the structures will be selected. Step indices can be
            negative, where negative index means counting reversely. E.g. "--save-as-poscars -2 -1 1 2 3" means saving
            the last two and first three steps.
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
- Output the trajectory of relaxation or MD
- Save the viberation modes

# How to build

- Make sure the public network is accesible to your machine
- If you still haven't installed Rust on your machine, follow [this](https://www.rust-lang.org/tools/install) to install
- Clone this repo, and `cd rsgrad`
- To run all tests: `cargo test`
- To get the binay: `cargo build --release`, and the `rsgrad` is on `rsgrad/target/release/rsgrad`
- Have fun!
