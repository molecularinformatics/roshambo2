# taken from https://github.com/pybind/cmake_example

import os
import re
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


# A CMakeExtension needs a sourcedir instead of a file list. Several may share
# one sourcedir; CMakeBuild.build_extensions() configures and builds it once.
class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


def _pybind11_cmake_dir() -> str:
    """
    Return the directory containing pybind11's CMake config files.
    Works inside PEP 517 build isolation.
    """
    try:
        import pybind11
        return pybind11.get_cmake_dir()
    except Exception:
        out = subprocess.check_output(
            [sys.executable, "-m", "pybind11", "--cmakedir"]
        )
        return out.decode().strip()


class CMakeBuild(build_ext):
    def build_extensions(self) -> None:
        # Both extensions come from the same CMakeLists.txt, and CMake already
        # builds every target in one configure+build pass. The stock
        # cmake_example pattern calls build_extension() once per Extension,
        # which ran the whole project twice here -- each pass compiling both
        # targets -- roughly doubling the build time.
        self.check_extensions_list(self.extensions)
        exts = list(self.extensions)
        if not exts:
            return
        sourcedir = exts[0].sourcedir
        assert all(e.sourcedir == sourcedir for e in exts), (
            "all CMakeExtensions must share the same sourcedir"
        )

        # Must be in this form due to bug in .resolve() only fixed in Python 3.10+
        # Both targets install to the same directory, so any one extension's
        # output path gives the right CMAKE_LIBRARY_OUTPUT_DIRECTORY.
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(exts[0].name)
        extdir = ext_fullpath.parent.resolve()

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            # Pass both the legacy (FindPythonInterp) and modern (FindPython)
            # hints. pybind11 >= 3 defaults to CMake's new FindPython, which
            # ignores the lowercase PYTHON_EXECUTABLE and would otherwise pick
            # a system interpreter whose dev headers may be absent -- which is
            # what happens under pip/uv build isolation.
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            f"-DEXAMPLE_VERSION_INFO={self.distribution.get_version()}",
            # Point CMake at pybind11's CMake config (fixes build isolation discovery)
            f"-Dpybind11_DIR={_pybind11_cmake_dir()}",
        ]
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja

                    ninja_executable_path = Path(ninja.BIN_DIR) / "ninja"
                    cmake_args += [
                        "-GNinja",
                        f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_executable_path}",
                    ]
                except ImportError:
                    pass

        else:
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        # One shared build tree for the whole CMake project (both targets).
        build_temp = Path(self.build_temp) / "roshambo2_cmake"
        build_temp.mkdir(parents=True, exist_ok=True)

        # Ensure CMake can also discover via prefix path
        env = os.environ.copy()
        p11 = _pybind11_cmake_dir()
        env["CMAKE_PREFIX_PATH"] = p11 + os.pathsep + env.get("CMAKE_PREFIX_PATH", "")

        subprocess.run(
            ["cmake", sourcedir, *cmake_args], cwd=build_temp, check=True, env=env
        )
        subprocess.run(
            ["cmake", "--build", ".", *build_args], cwd=build_temp, check=True, env=env
        )


setup(
    ext_modules=[CMakeExtension("_roshambo2_cuda"), CMakeExtension("_roshambo2_cpp")],
    cmdclass={"build_ext": CMakeBuild},
    scripts=[
        "roshambo2/scripts/prepare_dataset_from_sdf.py",
        "roshambo2/scripts/roshambo2_server_app.py",
    ],
)
