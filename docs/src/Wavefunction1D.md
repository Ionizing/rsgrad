# Wavefunction 1D

This command is similar to `rsgrad wav3d`. The only difference is that `wav1d` integrates
the wavefunction over specified plane, and produces the averaged one dimensional wavefunction.

This command is useful when you are looking for image potential states.

{{#include ./Wave1D-example-inline0.html}}

## Help Message
```shell
$ rsgrad wav1d --help
rsgrad-wav1d 
Plot wavefunction in realspace, then integrate over some plane, and save it as '.txt' file

USAGE:
    rsgrad wav1d [OPTIONS]

OPTIONS:
        --axis <AXIS>
            Integration direction. e.g. if 'z' is provided, the XoY plane is integrated
            
            [default: z]
            [possible values: X, Y, Z]

    -b, --ibands <IBANDS>...
            Select band index, starting from 1

    -d, --detail
            Show the eigen values and band occupations of current WAVECAR.
            
            This flag should be used with `--list`

        --gamma-half <GAMMA_HALF>
            Gamma Half direction of WAVECAR. You need to set this to 'x' or 'z' when processing
            WAVECAR produced by `vasp_gam`
            
            [possible values: x, z]

    -h, --help
            Print help information

        --htmlout <HTMLOUT>
            Specify the file name to be written with html wav1d data
            
            [default: wav1d.html]

    -k, --ikpoints <IKPOINTS>...
            Select kpoint index, starting from 1
            
            [default: 1]

    -l, --list
            List the brief info of current WAVECAR

    -s, --ispins <ISPINS>...
            Select spin index, starting from 1
            
            [default: 1]

        --scale <SCALE>
            Scale the wavefunction
            
            [default: 10]

        --show
            Open the browser and show the plot immediately

        --to-inline-html
            Render the plot and print thw rendered code to stdout

        --txtout <TXTOUT>
            Specify the file name to be written with raw wav1d data
            
            [default: wav1d.txt]

    -w, --wavecar <WAVECAR>
            WAVECAR file name
            
            [default: ./WAVECAR]
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
