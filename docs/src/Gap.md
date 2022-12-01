# Gap

This command (`rsgrad gap`) calculates the band gap and prints positions of VBM and CBM.

There are two ways to run this command:
1. Use `WAVECAR` only. All required information can be extracted from `WAVECAR`.
2. Or use `PROCAR` and `OUTCAR` instead. `PROCAR` is read for the band energies, band occupations and k-point coordinates;
`OUTCAR` is read for the Fermi level only.

Note: Usually it's faster to parse results from `WAVECAR`.

If your system is a non-SCF calculation for band structure, you may need to specify the correct Fermi level by
`rsgrad gap -e <put the correct E-fermi here>`, the correct Fermi level can be obtained from the OUTCAR of SCF calculation.


## Help Message
```shell
$ rsgrad gap --help
rsgrad-gap
Find band gap and print positions of VBM and CBM

USAGE:
    rsgrad gap [OPTIONS]

OPTIONS:
    -e, --efermi <EFERMI>      Specify Fermi level, if left empty, this value would be read from
                               WAVECAR or OUTCAR
    -h, --help                 Print help information
    -o, --outcar <OUTCAR>      OUTCAR file name, this file is parsed to get Fermi level only
                               [default: OUTCAR]
    -p, --procar <PROCAR>      PROCAR file name, OUTCAR is also needed to get Fermi level [default:
                               PROCAR]
    -w, --wavecar <WAVECAR>    WAVECAR file name, no more files needed [default: WAVECAR]
```

## Example

1. For metal system:
```shell
[2022-12-01T11:30:49Z INFO  rsgrad::commands::gap] Trying to parse "WAVECAR" ...
----------------------------------------
 Current system is         Metal
----------------------------------------
[2022-12-01T11:30:49Z INFO  rsgrad] Time used: 166.251241ms
```

2. For non-spin-polarized direct gap system:
```shell
[2022-12-01T11:37:15Z INFO  rsgrad::commands::gap] Trying to parse "WAVECAR" ...
--------------------------------------------------------------------------------
 Current system has    Direct Gap    of   0.328    eV
  CBM @ k-point     1 of ( 0.000, 0.000, 0.000) , band   819 of    0.322 eV
  VBM @ k-point     1 of ( 0.000, 0.000, 0.000) , band   818 of   -0.006 eV
--------------------------------------------------------------------------------
[2022-12-01T11:37:15Z INFO  rsgrad] Time used: 253.556682ms
```

3. For spin-polarized direct gap system:
```shell
$ rsgrad gap
[2022-12-01T11:38:16Z INFO  rsgrad::commands::gap] Trying to parse "WAVECAR" ...
--------------------------------------------------------------------------------
    ====================  Gap Info For     SPIN UP      ====================
 Current system has    Direct Gap    of   0.976    eV
  CBM @ k-point     1 of ( 0.000, 0.000, 0.000) , band   129 of    0.498 eV
  VBM @ k-point     1 of ( 0.000, 0.000, 0.000) , band   128 of   -0.478 eV
    ====================  Gap Info For    SPIN DOWN     ====================
 Current system has    Direct Gap    of   0.980    eV
  CBM @ k-point     1 of ( 0.000, 0.000, 0.000) , band   129 of    0.499 eV
  VBM @ k-point     1 of ( 0.000, 0.000, 0.000) , band   128 of   -0.481 eV
--------------------------------------------------------------------------------
[2022-12-01T11:38:16Z INFO  rsgrad] Time used: 31.744371ms
```

4. For non-spin-polarized indirect gap system:
```shell
$ rsgrad gap
[2022-12-01T11:40:52Z INFO  rsgrad::commands::gap] Trying to parse "WAVECAR" ...
--------------------------------------------------------------------------------
 Current system has   Indirect Gap   of   1.160    eV
  CBM @ k-point    25 of ( 0.000, 0.444, 0.444) , band     5 of    1.097 eV
  VBM @ k-point    17 of ( 0.000, 0.000, 0.000) , band     4 of   -0.062 eV
--------------------------------------------------------------------------------
[2022-12-01T11:40:52Z INFO  rsgrad] Time used: 22.565948ms
```
