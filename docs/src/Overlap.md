# Overlap

`rsgrad overlap` computes wavefunction overlaps between two `WAVECAR` files.

By default it includes the PAW one-centre correction and evaluates the all-electron overlap

`S_{ij}(k) = <Φ_i^a(k)|Φ_j^b(k)>`.

If you only want the pseudo-wavefunction overlap, use `--pseudo`.

## Help Message

```shell
$ rsgrad overlap --help
Compute AE wavefunction overlaps between two WAVECARs.

Computes S_{ij}(k) = <Φ_i^a(k)|Φ_j^b(k)> including the PAW one-centre correction.

Usage: rsgrad overlap [OPTIONS] --wavecar2 <WAVECAR2> --poscar2 <POSCAR2> --ibands <IBANDS>... --jbands <JBANDS>...

Options:
      --wavecar1 <WAVECAR1>
          First WAVECAR file path
          
          [default: WAVECAR]

      --wavecar2 <WAVECAR2>
          Second WAVECAR file path (required)

      --poscar1 <POSCAR1>
          First POSCAR file path
          
          [default: POSCAR]

      --poscar2 <POSCAR2>
          Second POSCAR file path (required)

      --potcar <POTCAR>
          POTCAR file path (shared between both structures)
          
          [default: POTCAR]

      --ispins <ISPINS>...
          Spin indices (1-indexed). Multiple values allowed
          
          [default: 1]

      --kpoints <KPOINTS>...
          K-point indices (1-indexed). Multiple values allowed
          
          [default: 1]

      --ibands <IBANDS>...
          Initial band index specs (e.g. "1:10" "15" "1:20:2"). 1-indexed

      --jbands <JBANDS>...
          Final band index specs (e.g. "1:10" "15" "1:20:2"). 1-indexed

      --pseudo
          If set, skip the PAW one-centre correction (pseudo-only overlap)

  -o, --output <OUTPUT>
          Output file path. Extension determines format: .npz → numpy, else text table
          
          [default: overlap.dat]

  -h, --help
          Print help (see a summary with '-h')
```

## Examples

Compute overlaps between two calculations and write a text table:

```shell
$ rsgrad overlap \
    --wavecar1 WAVECAR \
    --wavecar2 ../other/WAVECAR \
    --poscar1 POSCAR \
    --poscar2 ../other/POSCAR \
    --potcar POTCAR \
    --ibands 1:10 \
    --jbands 1:10 \
    -o overlap.dat
```

Write the overlap data as a NumPy archive:

```shell
$ rsgrad overlap \
    --wavecar2 ../other/WAVECAR \
    --poscar2 ../other/POSCAR \
    --ibands 1:20 \
    --jbands 1:20 \
    -o overlap.npz
```

Use `--pseudo` if you want pseudo-only overlaps without the PAW correction:

```shell
$ rsgrad overlap --wavecar2 ../other/WAVECAR --poscar2 ../other/POSCAR --ibands 1:10 --jbands 1:10 --pseudo
```
