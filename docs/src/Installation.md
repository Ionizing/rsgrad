# Installation

There are two ways to install `rsgrad`, downloading the pre-built binaries and building from scratch.

## Download the Binary from GitHub Release

We provide the pre-built binaries at [GitHub Release](https://github.com/Ionizing/rsgrad/releases/latest).

There are four platforms we support with pre-built binaries

- `rsgrad-<VERSION>-linux-x86_64-musl.tar.gz`: Statically linked binary. Can run on almost all Linux distros.
However, the `musl libc` may have performance issues in some scenes like `wave1d`, `wav3d` etc.

- `rsgrad-<VERSION>-linux-x86_64.tar.gz`: Dynamically linked against `glibc` and other fundamental
system libraries.  
For the binaries with _0.3.6_ci-test_ or upper version, they can run on most Linux distros
(e.g. CentOS 7, Ubuntu 18.04).  
If it still complains with message like `/lib64/libc.so.6: version 'GLIBC_2.18' not found`, just switch to `-musl`
binary or build from scratch.

- `rsgrad-<VERSION>-macos-x86_64.tar.gz`: Binary for macOS (Intel). `rsgrad` is developed on macOS 13, which
means you are don't need to worry about the availability if you are using macOS.

- `rsgrad-<VERSION>-macos-aarch64.tar.gz`: Binary for macOS (Apple Silicon).

- `rsgrad-<VERSION>-windows-x86_64.zip`: Binary for Windows system. We have no plan to support 32-bit systems.

Download the corresponding binary and put it in `$PATH` (or `%PATH%` on Windows), and have fun!

## Build from Scratch

__The build requires Internet accessibility.__

- If you haven't installed Rust toolchain, run `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`
and follow the instructions to configure it. When you finished the toolchain installation, you should be able to
run `cargo` in command-line. Detailed instructions can be found [here](https://rustup.rs/).

- If you have installed Rust toolchain already, just run `cargo install --git https://github.com/Ionizing/rsgrad`.
Several minutes later, the `rsgrad` should be installed to `~/.cargo/bin`, no need to modify the `$PATH` (or `%PATH%` on Windows).
