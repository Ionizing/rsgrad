# Model Non-adiabatic Coupling

This command (`rsgrad model-nac`) calculates the model non-adiabatic coupling
for NAMD-LMI (a subset of Hefei-NAMD).

**NOTE**: This command calculates momentum matrix element `<i|p|j>` using velocity gauge only,
thus the transition dipole moment in length gauge `<i|r|j>` is not used.


## Help Message

```shell
$ rsgrad model-nac --help
Calculate model non-adiabatic coupling (NAC) for NAMD-LMI (a subset of Hefei-NAMD).

This NAC contains eigenvalues and transition dipole moment (TDM) calculated from selected WAVECAR. The phonon contribution is wipped out, which means the eigenvalues and TDM stay staic over the time, and the NAC (<i|d/dt|j>) vanishes.

Detailed fields of the produced file:
 - ikpoint: K point index, counts from 1;
 - nspin: number of spin channels;
 - nbands: Number of total bands in WAVECAR;
 - brange: selected band indices, <start> ..= <end>, index counts from 1;
 - nbrange: number of selected bands, nbrange = end - start + 1;
 - efermi: Fermi's level;
 - potim: ionic time step;
 - temperature: 1E-6 Kelvin as default;
 - eigs: band eigenvalues;
 - pij_r/pij_i: real and imaginary part of <i|p|j>;
 - proj: Projection on each orbitals of selected bands, cropped from PROCAR.

Some fields not listed here are not meaningful but essential for the NAMD-LMI.

Usage: rsgrad model-nac [OPTIONS]

Options:
  -w, --wavecar <WAVECAR>
          WAVECAR file name

          [default: ./WAVECAR]

      --gamma-half <GAMMA_HALF>
          Gamma Half direction of WAVECAR. You need to set this to 'x' or 'z' when processing WAVECAR produced by `vasp_gam`

          [possible values: x, z]

  -k, --ikpoint <IKPOINT>
          One selected K-Point index, count starts from 1.

          Example: --ikpoint 2

          [default: 1]

      --brange <BRANGE> <BRANGE>
          Selected band range, starts from 1.

          Example: --brange 24 42

      --potim <POTIM>
          Ionic time step in femtosecond (fs)

          [default: 1]

  -p, --procar <PROCAR>
          PROCAR file name for band projection parsing

          [default: ./PROCAR]

  -o, --h5out <H5OUT>
          Output file name

          [default: ./NAC-0K.h5]

  -h, --help
          Print help (see a summary with '-h')
```

## Example

The `--brange <start> <end>` is required, or it will complain with error.

```shell
rsgrad model-nac --brange 107 110
[2024-11-18T04:17:11Z INFO  rsgrad::commands::modelnac] Reading WAVECAR: "./WAVECAR"
[2024-11-18T04:17:11Z INFO  rsgrad::commands::modelnac] Reading PROCAR: "./PROCAR"
[2024-11-18T04:17:11Z INFO  rsgrad::commands::modelnac] Saving to "./NAC-0K.h5"
[2024-11-18T04:17:11Z INFO  rsgrad] Time used: 87.687584ms
```
