# Installation

This software is distributed in source code form. It can be built with pip
inside an appropriate conda environment, or with
[uv](https://docs.astral.sh/uv/). Please follow the steps below.

## Installing with conda

1. Clone this repo (or extract from a tar archive/zip folder):
```
git clone <repo>
```

2. Create a conda environment with the provided yaml file:
```
cd roshambo2
conda env create -n roshambo2 -f environment.yaml
conda activate roshambo2
```

3. Install this package:
```
pip install .
```

## Installing with uv

[uv](https://docs.astral.sh/uv/) is a fast Python package and project manager.
If you have uv installed, you can install all dependencies and build the package
with a single command from the repo root:

```bash
git clone <repo>
cd roshambo2
uv sync
```

This reads `pyproject.toml`, creates a virtual environment, resolves and
installs all dependencies, and builds the C++/CUDA extensions. Activate
the environment with:

```bash
source .venv/bin/activate
```

## Building distributable wheels with uv

`uv sync` builds the extensions in place for development. To produce a standalone
wheel you can install elsewhere, use `uv build`. Build a wheel for the current
Python:

```bash
uv build --wheel
```

The wheel is written to `dist/`, e.g.
`roshambo2-1.0.0-cp311-cp311-linux_x86_64.whl`. Because the package contains
compiled C++/CUDA extensions, a wheel is specific to the Python version, OS, and
CPU architecture it was built on, and embeds code for the CUDA architecture(s)
selected at build time. `nvcc` must be on your `PATH` (see *Troubleshooting*).

> **Clean the build tree if a build fails.** CMake caches state under `build/`; a
> previous (especially failed) configure can poison the next one with errors like
> `generator : Ninja` or stale compiler/interpreter paths. If a build fails
> unexpectedly, run `rm -rf build dist` and retry.

### Alternate Python versions

Pass `-p`/`--python` with the desired version (uv fetches it if necessary). Each
build emits a wheel tagged for that interpreter (`cp311`, `cp312`, ...):

```bash
uv python install 3.11 3.12        # once, if those versions aren't available
uv build --wheel -p 3.11
uv build --wheel -p 3.12
```

Successful builds for different versions use version-specific subdirectories of
`build/` and do not collide, but `rm -rf build` between builds is the safe choice
if anything goes wrong.

### Different CUDA architectures

Select the target GPU architecture(s) with `CMAKE_CUDA_ARCHITECTURES` passed via
the `CMAKE_ARGS` environment variable (architecture values are in the table under
*CUDA build options* below). The default already targets sm_80 through sm_120
plus compute_120 PTX, so a default wheel is portable across current GPUs; narrow
it if you want a smaller wheel and know where it will run.

```bash
# a single architecture (e.g. H100)
CMAKE_ARGS="-DCMAKE_CUDA_ARCHITECTURES=90" uv build --wheel

# several, semicolon-separated, for a portable "fat" wheel
CMAKE_ARGS="-DCMAKE_CUDA_ARCHITECTURES=80;89;90;120" uv build --wheel
```

`CMAKE_ARGS` and `-p` combine, e.g. an sm_90 wheel for Python 3.12:

```bash
CMAKE_ARGS="-DCMAKE_CUDA_ARCHITECTURES=90" uv build --wheel -p 3.12
```

### Installing a built wheel

```bash
uv pip install dist/roshambo2-1.0.0-cp311-cp311-linux_x86_64.whl
# or into a specific interpreter / environment:
uv pip install --python /path/to/python dist/roshambo2-*.whl
```

### Building against an existing environment

For faster iteration you can build against the currently-active environment
(reusing its installed build tools and your system `cmake`) with
`--no-build-isolation`:

```bash
uv build --wheel --no-build-isolation
```

This requires `setuptools`, `pybind11`, `cmake` 3.24 or newer, and `ninja` to be
installed in that environment.

## Documentation install

To build and view the documentation please go into the `doc` folder.
You will then need to install the extra dependencies (these are listed in doc/env.yaml) to your conda environment:

```
conda install  -c sphinx pydata-sphinx-theme myst-parser
```
You can then build the html docs with:
```
make html
```
And you can view them in your browser:
```
xdg-open build/html/index.html 
```


## Requirements

A CMake 3.24 or later build system is used to compile the C++ and CUDA code.
You will need a working C++ compiler, `nvcc`, and CMake 3.24+.
CMake will need to be able to find your CUDA installation.

## CUDA build options

### Target GPU architecture

By default the build targets `80-real;86-real;89-real;90-real;120-real;120-virtual`:
native SASS for sm_80, sm_86, sm_89, sm_90 and sm_120, plus compute_120 PTX so
the driver can JIT for any GPU newer than those. The result runs on any current
GPU without a rebuild, at the cost of compile time and binary size.

To target a different set of architectures — for example to cut build time, or
to reach a GPU older than sm_80 — pass `CMAKE_CUDA_ARCHITECTURES` via the
`CMAKE_ARGS` environment variable:

```bash
# Single arch
CMAKE_ARGS="-DCMAKE_CUDA_ARCHITECTURES=90" pip install .

# Multiple arches (semicolon-separated)
CMAKE_ARGS="-DCMAKE_CUDA_ARCHITECTURES=80;89;90" pip install .
```

The same flag works with `uv`:

```bash
CMAKE_ARGS="-DCMAKE_CUDA_ARCHITECTURES=80;89;90" uv sync
```

For fast local iteration, `-DROSHAMBO_ARCH_NATIVE=ON` compiles only for the
GPU(s) in the build machine, which is much quicker:

```bash
CMAKE_ARGS="-DROSHAMBO_ARCH_NATIVE=ON" pip install .
```

Do not use it for anything you intend to copy elsewhere: the extension will not
load on a different architecture.

Common architecture values by GPU, with corresponding cloud instance types:

| GPU  | Architecture | `CMAKE_CUDA_ARCHITECTURES` | AWS instance | GCP instance |
|------|-------------|:--------------------------:|--------------|--------------|
| P100                        | Pascal       | `60`  | —         | N1 |
| P4                          | Pascal       | `61`  | —         | N1 |
| V100                        | Volta        | `70`  | —         | N1 |
| T4 / T4g                    | Turing       | `75`  | G4dn, G5g | T4 |
| A100                        | Ampere       | `80`  | P4        | A2 |
| A10g                        | Ampere       | `86`  | G5        | — |
| L4/L40S                     | Ada Lovelace | `89`  | G6, G6e   | G2 |
| H100 / H200                 | Hopper       | `90`  | P5        | A3 |
| B200 / GB200 / GB300        | Blackwell    | `100` | P6, P6e   | A4 |
| RTX PRO 4500/6000 Blackwell | Blackwell    | `120` | G7, G7e   | G4 |

Architectures below `80` are not in the default list and must be requested
explicitly. Note that CUDA 13 dropped support for `compute_60` altogether, so
`60` only works with an older toolkit.

If you are unsure of your GPU's compute capability, run:
```bash
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader
```

## Troubleshooting

If cmake complains about not finding cuda you will need to make sure `nvcc` can be found on your path.
You might need to do something like this:
```bash
export PATH=/usr/local/cuda/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH}
```
