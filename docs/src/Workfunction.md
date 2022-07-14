# Workfunction

`rsgrad workfunc` reads _LOCPOT_ to calculate the work function. For further discussion about
work function, please read [this blog](https://ionizing.page/post/vasp-dipol-correction-work-function/).

{{#include ./Workfunction-example0.html}}

__Note: `rsgrad workfunc` does not show the vacuum level directly, hence you need to find it by
hovering on the step of the work function. Take the work function above as an example, the
vacuum level should be about 4.048eV around 40Å.__

## Help Message

```shell
$ rsgrad workfunc --help
rsgrad-workfunc
Calculate work-function from LOCPOT file, OUTCAR is also needed to get the Fermi level.

The work function is calculated by plannar integration of the data cube in LOCPOT. The selected axis
should be perpendicular to the other two axises.

USAGE:
    rsgrad workfunc [OPTIONS] [LOCPOT]

ARGS:
    <LOCPOT>
            LOCPOT file path. Turn on 'LVHAR' in INCAR to get the electro-static potential saved it

            [default: ./LOCPOT]

OPTIONS:
        --axis <AXIS>
            Integration direction. e.g. if 'z' is provided, the XoY plane is integrated

            [default: z]
            [possible values: X, Y, Z]

    -h, --help
            Print help information

    -o, --htmlout <HTMLOUT>
            Write the plot to html and view it in the web browser

            [default: ./workfunction.html]

        --outcar <OUTCAR>
            OUTCAR file path. This file is needed to get the E-fermi level and lattice properties

            [default: ./OUTCAR]

        --show
            Open default browser to see the plot immediately

        --to-inline-html
            Render the plot and print the rendered code to stdout

        --txtout <TXTOUT>
            Write the raw plot data as txt file in order to replot it with more advanced tools

            [default: locpot.txt]
```

## Example

In general, you don't need t provide additional arguments.

```shell
$ rsgrad workfunc
[2022-07-14T14:41:38Z INFO  rsgrad::commands::workfunc] Reading "./OUTCAR"
[2022-07-14T14:41:38Z INFO  rsgrad::commands::workfunc] Reading electro-static potential data from "./LOCPOT"
[2022-07-14T14:41:39Z INFO  rsgrad::commands::workfunc] Writing raw plot data to "locpot.txt"
[2022-07-14T14:41:39Z INFO  rsgrad::commands::workfunc] Writing to "./workfunction.html"
[2022-07-14T14:41:39Z INFO  rsgrad] Time used: 1.352347116s
```

{{#include ./Workfunction-example1.html}}
