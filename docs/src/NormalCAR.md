# NormalCAR

`rsgrad normalcar` writes PAW projector coefficients from `WAVECAR`, `POSCAR`, and `POTCAR`.

The generated `NormalCAR` is useful for workflows that need the projector coefficients
`β = ⟨p̃_i|ψ̃⟩`, such as spin-orbit and other PAW-based post-processing.

## Help Message

```shell
$ rsgrad normalcar --help
Write NormalCAR (PAW projector coefficients) from POSCAR + POTCAR + WAVECAR.

The NormalCAR stores the PAW projector coefficients β = ⟨p̃_i|ψ̃⟩ for every spin/k-point/band, useful for spinorb and other post-processing tools.

Usage: rsgrad normalcar [OPTIONS]

Options:
      --wavecar <WAVECAR>
          WAVECAR file path
          
          [default: WAVECAR]

      --poscar <POSCAR>
          POSCAR file path
          
          [default: POSCAR]

      --potcar <POTCAR>
          POTCAR file path
          
          [default: POTCAR]

  -o, --output <OUTPUT>
          Output file path. If extension is .npz, write numpy npz format; otherwise write binary NormalCAR
          
          [default: NormalCAR]

  -h, --help
          Print help (see a summary with '-h')
```

## Examples

Write the binary `NormalCAR` used by VASP-style post-processing:

```shell
$ rsgrad normalcar --wavecar WAVECAR --poscar POSCAR --potcar POTCAR -o NormalCAR
```

Write the same coefficients as a NumPy archive:

```shell
$ rsgrad normalcar --wavecar WAVECAR --poscar POSCAR --potcar POTCAR -o cproj.npz
```

If the output file ends with `.npz`, rsgrad writes a NumPy-friendly archive. Otherwise it writes
the binary `NormalCAR` format.
