# Transition Dipole Moment

This command (`rsgrad tdm`) calculates the transition dipole moment(tdm)
between given bands, and then plot the peaks with plotly.

{{#include ./tdm-inline.html}}

Note: This command can only calculate the TDM between bands in same k-point. The inter-kpoint
transition is not supported yet. Also, this commands calculates the TDM in reciprocal space by

\\[
TDM_{i\to j} = \langle \phi_j | e\cdot\vec{r} | \phi_i \rangle = e\cdot\frac{i\hbar}{m_e\Delta E_{ij}} \langle \phi_j | \vec{p} | \phi_i \rangle
\\]

Also, this command only calculates the TDM if \\(j > i\\), otherwise it will ignore the invalid \\(i\to j\\) pairs.

## Help Message
```shell
$ rsgrad tdm --help
rsgrad-tdm
Calculate Transition Dipole Moment (TDM) between given bands.

Note: This command can only calculate the TDM between bands in same k-point. The inter-kpoint
transition is not supported yet. Also, this commands calculates the TDM in reciprocal space by

tdm_{i->j} = <phi_j|e*r|phi_i> = i*ħ/(ΔE*m)*<phi_j|p|phi_i>

USAGE:
    rsgrad tdm [OPTIONS] --ibands <IBANDS>... --jbands <JBANDS>...

OPTIONS:
        --barwidth <BARWIDTH>
            Specify the width of bars in the center of peaks. (eV)

            [default: 0.1]

        --gamma-half <GAMMA_HALF>
            Gamma Half direction of WAVECAR. You need to set this to 'x' or 'z' when processing
            WAVECAR produced by `vasp_gam`

            [possible values: x, z]

    -h, --help
            Print help information

        --htmlout <HTMLOUT>
            Write the plot of TDM to html file

            [default: tdm_smeared.html]

    -i, --ibands <IBANDS>...
            Initial band indices, start from 1

    -j, --jbands <JBANDS>...
            Final band indices, starts from 1

    -k, --ikpoint <IKPOINT>
            K-point index, starts from 1

            [default: 1]

        --peakout <PEAKOUT>
            Write the TDM peaks to raw txt file

            [default: tdm_peaks.txt]

    -s, --ispin <ISPIN>
            Spin index, 1 for up, 2 for down

            [default: 1]
            [possible values: 1, 2]

        --show
            Open the default browser to show the plot

        --sigma <SIGMA>
            Smearing width, in eV

            [default: 0.05]

        --to-inline-html
            Print the inline HTML to stdout

        --txtout <TXTOUT>
            Write the summed and smeared TDM to raw txt file

            [default: tdm_smeared.txt]

    -v, --verbose
            Print the calculated TDM to screen

    -w, --wavecar <WAVECAR>
            WAVECAR file path

            [default: ./WAVECAR]
```

## Example

The `-i` and `-j` arguments are required, or it will complain with error.

```shell
$ rsgrad tdm -i 2 -j {3..5}
[2022-08-10T07:26:46Z INFO  rsgrad::commands::tdm] Reading "./WAVECAR"
[2022-08-10T07:26:46Z INFO  rsgrad::commands::tdm] Writing peaks info to "tdm_peaks.txt" ...
[2022-08-10T07:26:46Z INFO  rsgrad::commands::tdm] Writing smeared TDM data to "tdm_smeared.txt"
[2022-08-10T07:26:46Z INFO  rsgrad::commands::tdm] Writing plot to "tdm_smeared.html"
[2022-08-10T07:26:46Z INFO  rsgrad] Time used: 56.655387ms
```

It produces the HTML like:

{{#include ./tdm-inline2.html}}

You can view the \\(x\\), \\(y\\) and \\(z\\) components of each peak by hovering on the thin bars
at the center of peaks.

**Note: If two peaks are in the same x position, one of them will be covered.**

Then the `tdm_peaks.txt` should be like:

```
# iband jband     E_i     E_j      ΔE        Tx       Ty       Tz
      2     3  -5.974  -3.675   2.298     1.311    2.556    0.000
      2     4  -5.974  -3.675   2.298     2.556    1.311    0.000
      2     5  -5.974  -0.975   4.999     0.000    0.000    5.309
```

And the `tdm_smeared.txt` should be like:

```
# E(eV)    Tx(Debye)   Ty(Debye)   Tz(Debye)
         0.298325         0.007692         0.007692         0.001912         0.017295
         0.300240         0.007706         0.007706         0.001914         0.017326
         0.302155         0.007721         0.007721         0.001915         0.017357
         0.304070         0.007736         0.007736         0.001917         0.017389
         0.305985         0.007751         0.007751         0.001918         0.017420
         0.307900         0.007766         0.007766         0.001920         0.017451
...
...
```

Then you can plot the TDM peaks with other tools like OriginPro.

If your `WAVECAR` is produced by `vasp_gam`, you'd better specify the gamma half direction by `--gamma-half x` or `--gamma-half z`,
or `rsgrad tdm` will produce wrong results.
