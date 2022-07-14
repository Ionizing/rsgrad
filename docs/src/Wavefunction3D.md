# Wavefunction 3D

Main functions of `rsgrad wav3d`
- List the brief information of _WAVECAR_;
- Save the selected wavefunction to _.vasp_ file.

## Help Message

```shell
$ rsgrad wav3d -h
rsgrad-wav3d
Plot wavefunction in realspace, and save it as '.vasp' file

USAGE:
    rsgrad wav3d [OPTIONS]

OPTIONS:
    -b, --ibands <IBANDS>...
            Select band index, starting from 1

    -e, --show-eigs-suffix
            Add eigen value suffix to the filename

        --gamma-half <GAMMA_HALF>
            Gamma Half direction of WAVECAR. You need to set this to 'x' or 'z' when processing
            WAVECAR produced by `vasp_gam` [possible values: x, z]

    -h, --help
            Print help information

    -k, --ikpoints <IKPOINTS>...
            Select kpoint index, starting from 1 [default: 1]

    -l, --list
            List the brief info of current WAVECAR

    -o, --output-parts <OUTPUT_PARTS>...
            Specify output part of the wavefunction [possible values: normsquared, ns, real, re,
            imag, im, reim]

    -p, --poscar <POSCAR>
            POSCAR filename, POSCAR is needed to get the real-space wavefunction [default: ./POSCAR]

        --prefix <PREFIX>
            Prefix of output filename [default: wav]

    -s, --ispins <ISPINS>...
            Select spin index, starting from 1 [default: 1]

    -w, --wavecar <WAVECAR>
            WAVECAR file name [default: ./WAVECAR]
```

## Example

- List the brief information of _WAVECAR_
