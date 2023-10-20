# Pseudo Potential

This command generates _POTCAR_ according to _POSCAR_. The usage is simple: `rsgrad pot` in the directory with _POSCAR_ presents.

To use this command, a configuration file is required like:

```
[functional-path]
PAW_PBE = "<path of PAW_PBE>"
PAW_LDA = "<path of PAW_LDA>"
```

and you need to put it at `~/.rsgrad` for Linux/macOS,
and `C:\Users\<YourUserName>\.rsgrad.toml` for Windows.

**Note**: `rsgrad pot` will give some hint if you have no idea on how to write the file.

Then you can specify the element type in _POSCAR_, for example:
```
Cd2 S2 HSE optimized
   1.00000000000000
     4.1704589747753724    0.0000000049514629    0.0000000000000000
    -2.0852299788115007    3.6117232476351271   -0.0000000000000000
     0.0000000000000000   -0.0000000000000000    6.7822968954189937
   Cd   S
     2     2
Direct
  0.6666669999999968  0.3333330000000032  0.5001535716160643
  0.3333330000000032  0.6666669999999968  0.0001535716160644
  0.6666669999999968  0.3333330000000032  0.8768464283839381
  0.3333330000000032  0.6666669999999968  0.3768464283839381
```

Finally `rsgrad-pot` will read the line `"   Cd   S "`, and find the corresponding _POTCAR_ and merge them.

In some cases, several elements does not have "standard" _POTCAR_ that without suffix. E.g. the `K/POTCAR`
is not available for PBE functional. We provide two solutions for it:
- Write `K_sv` or other available K-related functionals in _POSCAR_
  ```
  Cd2 S2 HSE optimized
     1.00000000000000
       4.1704589747753724    0.0000000049514629    0.0000000000000000
      -2.0852299788115007    3.6117232476351271   -0.0000000000000000
       0.0000000000000000   -0.0000000000000000    6.7822968954189937
     Cd   K_sv
       2     2
  Direct
    0.6666669999999968  0.3333330000000032  0.5001535716160643
    ...
  ```
- Use `aliases` field in `.rsgrad.toml`
  ```
  [functional-path]
  ...
  aliases = { K = "K_sv" }
  ```
  or
  ```
  [functional-path]
  ...
  [functional-path.aliases]
  K = "K_sv"
  ```

## Help Message

```
rsgrad-pot
Generate the POTCAR according to POSCAR

USAGE:
    rsgrad pot [OPTIONS] [FUNCTIONAL]

ARGS:
    <FUNCTIONAL>
            Specify the functional type, now only "PAW_PBE"(or "paw_pbe") and "PAW_LDA"(or
            "paw_lda") are available

            [default: PAW_PBE]

OPTIONS:
    -c, --config <CONFIG>
            Specify the configuration file, if left blank, rsgrad will read `.rsgrad.toml` at your
            home dir

    -h, --help
            Print help information

    -p, --poscar <POSCAR>
            Specify the POSCAR file

            Note: the elements symbols' line should involve the specified valence configuration if
            you need. For example, if you need a `K_sv` POTCAR for `K`, just write `K_sv` in the
            POSCAR

            [default: ./POSCAR]

    -s, --save-in <SAVE_IN>
            Specify where the `POTCAR` would be written

            [default: ./]
```

## Example

```shell
$ rsgrad pot
[2023-10-20T16:19:48Z INFO  rsgrad::settings] Reading rsgrad settings from "/Users/ionizing/.rsgrad.toml" ...
[2023-10-20T16:19:48Z INFO  rsgrad::settings] Checking rsgrad setting availability ...
[2023-10-20T16:19:48Z INFO  rsgrad::commands::pot] Reading POSCAR file "../tests/POSCAR.CdS_HSE" ...
[2023-10-20T16:19:48Z INFO  rsgrad::vasp_parsers::potcar] Reading POTCAR from "/Users/ionizing/apps/vasp/pp/potpaw_PBE.54/Cd/POTCAR"
[2023-10-20T16:19:48Z INFO  rsgrad::vasp_parsers::potcar] Reading POTCAR from "/Users/ionizing/apps/vasp/pp/potpaw_PBE.54/S/POTCAR"
[2023-10-20T16:19:48Z INFO  rsgrad::commands::pot] Writing POTCAR to "./POTCAR" ...
[2023-10-20T16:19:48Z INFO  rsgrad] Time used: 5.524643ms
```
