# UnitConversion

This command (`rsgrad uc`) converts input quantity to units representing same energy.

Metric prefixes are supported, from _Atto_(a, 10<sup>-18</sup>) to _Exa_(E, 10<sup>18</sup>).

Here are supported prefixes and corresponding acronyms:

| Metric prefix   | Acronym   | Base 10          |
| --------------- | --------- | ---------------- |
| atto            | a         | 10<sup>-18</sup> |
| femto           | f         | 10<sup>-15</sup> |
| pico            | p         | 10<sup>-12</sup> |
| nano            | n         | 10<sup>-9</sup>  |
| micro           | μ         | 10<sup>-6</sup>  |
| milli           | m         | 10<sup>-3</sup>  |
| --              | --        | 10<sup>0</sup>   |
| Kilo            | K         | 10<sup>3</sup>   |
| Mega            | M         | 10<sup>6</sup>   |
| Giga            | G         | 10<sup>9</sup>   |
| Tera            | T         | 10<sup>12</sup>  |
| Peta            | P         | 10<sup>15</sup>  |
| exa             | E         | 10<sup>18</sup>  |

If you use acronyms, note that prefixes smaller than 1 (from _a_, atto to _m_, micro) must be
lowercased whereas the prefixes lager than 1 (from _K_, Kilo to _E_, Exa) must be uppercased.

Supported units:

| Unit                | Symbol  |
| ------------------- | -------:|
| Hartree             | Ha      |
| Wavenumber          | cm-1    |
| Temperature(Kelvin) | K       |
| Calorie per mole    | Cal/mol |
| Joule per mole      | J/mol   |
| Wavelength          | m       |
| Period              | s       |
| Electron volt       | eV      |
| Frequency           | Hz      |

## Help Message
```shell
$ rsgrad uc
Conversion between various energy units

Usage: rsgrad uc [INPUT]...

Arguments:
  [INPUT]...  Input energy quantity to be converted. Multiple input are supported

Options:
  -h, --help  Print help

Try `rsgrad uc 298K` to see what happens.
```

## Examples

```shell
$ rsgrad uc 298K
==================== Processing input "298K" ====================
  298.000000 K ==   943.709369 μHa
  298.000000 K ==   207.125149 cm-1
  298.000000 K ==   298.000000 K
  298.000000 K ==   592.202785 Cal/mol
  298.000000 K ==     2.477776 KJ/mol
  298.000000 K ==    48.281101 μm
  298.000000 K ==   161.048425 fs
  298.000000 K ==    25.679653 meV
  298.000000 K ==     6.209312 THz
================================================================================

[2024-09-11T07:21:07Z INFO  rsgrad] Time used: 470.398µs
```


**You can also input multiple quantities at a time**

```shell
$ rsgrad uc 298K 25meV 300THz
==================== Processing input "298K" ====================
  298.000000 K ==   161.048425 fs
  298.000000 K ==   207.125149 cm-1
  298.000000 K ==     2.477776 KJ/mol
  298.000000 K ==    25.679653 meV
  298.000000 K ==    48.281101 μm
  298.000000 K ==   943.709369 μHa
  298.000000 K ==   592.202785 Cal/mol
  298.000000 K ==     6.209312 THz
  298.000000 K ==   298.000000 K
================================================================================

==================== Processing input "25meV" ====================
   25.000000 meV ==   165.426708 fs
   25.000000 meV ==   201.643250 cm-1
   25.000000 meV ==     2.412198 KJ/mol
   25.000000 meV ==    25.000000 meV
   25.000000 meV ==    49.593677 μm
   25.000000 meV ==   918.732590 μHa
   25.000000 meV ==   576.529191 Cal/mol
   25.000000 meV ==     6.044973 THz
   25.000000 meV ==   290.112953 K
================================================================================

==================== Processing input "300THz" ====================
  300.000000 THz ==     3.333333 fs
  300.000000 THz ==    10.007154 Kcm-1
  300.000000 THz ==   119.712599 KJ/mol
  300.000000 THz ==     1.240700 eV
  300.000000 THz ==   999.308150 nm
  300.000000 THz ==    45.594872 mHa
  300.000000 THz ==    28.611998 KCal/mol
  300.000000 THz ==   300.000000 THz
  300.000000 THz ==    14.397729 KK
================================================================================

[2024-09-11T07:27:53Z INFO  rsgrad] Time used: 544.171µs
```
