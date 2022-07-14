# POTCAR

This command (`rsgrad pot`) is quite simple: it reads _POSCAR_ and generates corresponding POTCAR.

`rsgrad pot` need a configuration file defines the directories of pseudo potential files. By default,
 the configuration file is located at `~/.rsgrad.toml`, showing

```
$ cat ~/.rsgrad.toml
[functional-path]
PAW_PBE = "~/apps/pp/potpaw_PBE.54"
PAW_LDA = "~/apps/pp/potpaw_LDA.54"
```

However, if you prefer other file names or paths, you may need to specify the configuration file manually
every time:

```shell
rsgrad pot -c /path/to/.rsgrad.toml
```

__Note: This command does NOT depend on Linux/BSD system's built-in `cat` binary, which means it can also
run on Windows.__

## Help Message
```shell
$ rsgrad pot --help
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

The reference _POSCAR_:

```
C  O
 1.0000000
      20.000000000       0.000000000       0.000000000
       0.000000000      20.000000000       0.000000000
       0.000000000       0.000000000      20.000000000
    C_h    O_h
      1      2
Direct
      0.5000000000      0.5000000000     -0.5000000000 !    C_h-001    1
      0.5754584708      0.4974961803      0.5000000000 !    O_h-001    2
      0.4249089652      0.5078477016      0.5000000000 !    O_h-002    3
```

Then run the command

```shell
$ rsgrad pot -p POSCAR_new
[2022-07-14T16:43:39Z INFO  rsgrad::settings] Reading rsgrad settings from "/Users/ionizing/.rsgrad.toml" ...
[2022-07-14T16:43:39Z INFO  rsgrad::settings] Checking rsgrad setting availability ...
[2022-07-14T16:43:39Z INFO  rsgrad::commands::pot] Reading POSCAR file "POSCAR_new" ...
[2022-07-14T16:43:39Z INFO  rsgrad::vasp_parsers::potcar] Reading POTCAR from "/Users/ionizing/apps/pp/potpaw_PBE.54/C_h/POTCAR"
[2022-07-14T16:43:39Z INFO  rsgrad::vasp_parsers::potcar] Reading POTCAR from "/Users/ionizing/apps/pp/potpaw_PBE.54/O_h/POTCAR"
[2022-07-14T16:43:39Z INFO  rsgrad::commands::pot] Writing POTCAR to "./POTCAR" ...
[2022-07-14T16:43:39Z INFO  rsgrad] Time used: 8.556795ms
```

We can verify the _POTCAR_:

```shell
$ grep TITEL POTCAR
   TITEL  = PAW_PBE C_h 06Feb2004
   TITEL  = PAW_PBE O_h 06Feb2004
```

The `C_h` and `O_h` are exactly consistent with the _POSCAR_'s element labels.
