![Rust](https://github.com/Ionizing/rsgrad/workflows/Rust/badge.svg)

Tracking the relaxation or MD progress of VASP calculation.

# Usage

## Main interface

```
$ ./rsgrad
rsgrad 0.1.0
@Ionizing  https, //github.com/Ionizing/rsgrad
A tool used to tracking the relaxation or MD progress of VASP calculation

USAGE:
    rsgrad [input] <SUBCOMMAND>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

ARGS:
    <input>    Specify the input OUTCAR file name [default: ./OUTCAR]

SUBCOMMANDS:
    help    Prints this message or the help of the given subcommand(s)
    list    Lists the brief info of current OUTCAR
    rlx     Tracking info associated with relaxation stuff
    trj     Operations about relaxation/MD trajectory
    vib     Tracking info associated with vibration stuff```
```

## List the brief info of current OUTCAR

```
$ ./rsgrad some_OUTCAR list
[2021-02-24T14:53:29Z INFO  rsgrad] Parsing input file "OUTCAR_ispin2" ...
    IBRION =          1
     NKPTS =         41
     NIONS =          3
       NSW =          3
     ISPIN =          2
   LSORBIT =      false
    EFERMI =    -2.2691
    NBANDS =         16
[2021-02-24T14:53:29Z INFO  rsgrad] Time used: 69.130887ms
```

## See the relaxation info of OUTCAR

```
$ ./rsgrad some_OUTCAR rlx
[2021-02-24T14:54:13Z INFO  rsgrad] Parsing input file "OUTCAR_ispin2" ...
  #Step  TOTEN_z/eV LgdE   Fmax #SCF Time/m Mag/muB
      1   -18.95729  1.3  0.000   27   0.48   0.600
      2   -18.95789 -3.2  0.000    6   0.11   0.600
      3   -18.95797 -4.1  0.000    4   0.07   0.600
[2021-02-24T14:54:13Z INFO  rsgrad] Time used: 70.359006ms
```

Full functions of `rsgrad rlx`

```
$ ./rsgrad some_OUTCAR rlx -h
rsgrad-trj 0.1.0
Operations about relaxation/MD trajectory

USAGE:
    rsgrad trj [FLAGS] [OPTIONS]

FLAGS:
    -h, --help               Prints help information
    -p, --save-as-poscars    Saves structures of given steps as POSCARs
    -d, --save-as-xdatcar    Saves total trajectory in XDATCAR format
    -x, --save-as-xsfs       Saves structures of given steps as XSFs
    -V, --version            Prints version information

OPTIONS:
        --save-in <save-in>                     Defines where the files would be saved [default: .]
    -i, --select-indices <select-indices>...    Selects the indices to operate
```


## Operations on vibration modes

```
$ ./rsgrad some_OUTCAR vib -l
[2021-02-24T14:58:36Z INFO  rsgrad] Parsing input file "OUTCAR_vibrations" ...
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
```

Full usage

```
$ ./rsgrad some_OUTCAR vib -h
rsgrad-vib 0.1.0
Tracking info associated with vibration stuff

USAGE:
    rsgrad vib [FLAGS] [OPTIONS]

FLAGS:
    -h, --help            Prints help information
    -l, --list            Shows vibration modes in brief
    -x, --save-as-xsfs    Saves each selected modes to XSF file
    -V, --version         Prints version information

OPTIONS:
        --save-in <save-in>                     Define where the files would be saved [default: .]
    -i, --select-indices <select-indices>...    Selects the indices to operate
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
- [X] A prettier output layout
- [X] Output the trajectory of relaxation or MD
- [ ] Display the unconverged atoms (will be implemented in the near future)
- [X] Save the viberation modes
- [X] More detailed error messages

# How to build

- Make sure the public network is accesible to your machine
- If you still haven't installed Rust on your machine, follow [this](https://www.rust-lang.org/tools/install) to install
- Clone this repo, and `cd rsgrad`
- To run all tests: `cargo test`
- To get the binay: `cargo build --release`, and the `rsgrad` is on `rsgrad/target/release/rsgrad`
- Have fun!
