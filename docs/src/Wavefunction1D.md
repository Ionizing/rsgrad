# Wavefunction 1D

This command is similar to `rsgrad wav3d`. The only difference is that `wav1d` integrates
the wavefunction over specified plane, and produces the averaged one dimensional wavefunction.

It also supports all-electron reconstruction through `--ae`.

This command is useful when you are looking for image potential states.

{{#include ./Wave1D-example-inline0.html}}

## Help Message
```shell
$ rsgrad wav1d --help
Plot wavefunction in realspace, then integrate over some plane, and save it as '.txt' file

Usage: rsgrad wav1d [OPTIONS]

Options:
  -w, --wavecar <WAVECAR>
          WAVECAR file name
          
          [default: ./WAVECAR]

  -s, --ispins [<ISPINS>...]
          Select spin index, starting from 1
          
          [default: 1]

  -k, --ikpoints [<IKPOINTS>...]
          Select kpoint index, starting from 1.
          
          You can input range directly: `-k 1..5 8..10`
          
          [default: 1]

  -b, --ibands [<IBANDS>...]
          Select band index, starting from 1.
          
          You can input range directly: `-b 5..10 14..19`

  -l, --list
          List the brief info of current WAVECAR

  -d, --detail
          Show the eigen values and band occupations of current WAVECAR.
          
          This flag should be used with `--list`

      --gamma-half <GAMMA_HALF>
          Gamma Half direction of WAVECAR. You need to set this to 'x' or 'z' when processing WAVECAR produced by `vasp_gam`
          
          [possible values: x, z]

      --poscar <POSCAR>
          POSCAR filename, required when reconstructing all-electron wavefunctions
          
          [default: ./POSCAR]

      --potcar <POTCAR>
          POTCAR filename, required when reconstructing all-electron wavefunctions
          
          [default: ./POTCAR]

      --ae
          Reconstruct the all-electron wavefunction density instead of using the pseudo-wavefunction

      --aecut <AECUT>
          AE energy cutoff in eV. Negative values mean |aecut| * pscut
          
          [default: -2]

      --txtout <TXTOUT>
          Specify the file name to be written with raw wav1d data
          
          [default: wav1d.txt]

      --htmlout <HTMLOUT>
          Specify the file name to be written with html wav1d data
          
          [default: wav1d.html]

      --axis <AXIS>
          Integration direction. e.g. if 'z' is provided, the XoY plane is integrated
          
          [default: z]
          [possible values: x, y, z]

      --to-inline-html
          Render the plot and print thw rendered code to stdout

      --show
          Open the browser and show the plot immediately

      --scale <SCALE>
          Scale the wavefunction
          
          [default: 10]

  -h, --help
          Print help (see more with '--help')
```


## Example

The usage is similar to `rsgrad wav3d`, too.

```shell
$ rsgrad wav1d -b {20..27}
[2022-07-14T12:54:35Z INFO  rsgrad::commands::wav1d] Reading WAVECAR: "./WAVECAR"
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Processing spin 1, k-point   1, band   20 ...
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Processing spin 1, k-point   1, band   24 ...
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Processing spin 1, k-point   1, band   22 ...
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Processing spin 1, k-point   1, band   23 ...
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Processing spin 1, k-point   1, band   21 ...
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Processing spin 1, k-point   1, band   26 ...
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Processing spin 1, k-point   1, band   27 ...
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Processing spin 1, k-point   1, band   25 ...
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Writing to "wav1d.html"
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Writing to "wav1d.txt"
[2022-07-14T12:54:36Z INFO  rsgrad::commands::wav1d] Printing inline html to stdout ...
[2022-07-14T12:54:36Z INFO  rsgrad] Time used: 864.911523ms
```

Then produces

{{#include ./Wave1D-example-inline.html}}

__Note: The Fermi level is shiftted to 0.0 eV__

## All-Electron Reconstruction

`rsgrad aewfc` has been merged into `rsgrad wav1d` for 1D integrated output.

Use `--ae --poscar POSCAR --potcar POTCAR` to integrate the reconstructed all-electron density:

```shell
$ rsgrad wav1d --ae --poscar POSCAR --potcar POTCAR -b 20..27
```

Current limitation: `--ae` is not available for non-collinear WAVECAR.
