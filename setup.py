
# adapted from https://github.com/pybind/cmake_example/blob/master/setup.py
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

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        # Must be in this form due to bug in .resolve() only fixed in Python 3.10+
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)  # type: ignore[no-untyped-call]
        extdir = ext_fullpath.parent.resolve()

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        print(f"Running build in {cfg} mode")

        # this is needed for pip installed version of pybind11
        try:
            import pybind11
            pybind11_dir = "-Dpybind11_DIR={}".format(pybind11.get_cmake_dir())
        except ImportError as e:
            print("pybind11 is a required module")

        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            f"-DBUILD_PYCGMF=On",
            f"-DCMAKE_BUILD_TYPE=Release",
            f"-Dcgmf.x.MPI=On",
            pybind11_dir,
        ]
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # In this example, we pass in the version to C++. You might not need to.
        cmake_args += [f"-DOMPLIBPY_VERSION_INFO={self.distribution.get_version()}"]  # type: ignore[attr-defined]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja  # noqa: F401

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

        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        subprocess.run(
            ["cmake", ext.sourcedir] + cmake_args, cwd=build_temp, check=True
        )
        subprocess.run(
            ["cmake", "--build", "."] + build_args, cwd=build_temp, check=True
        )


def install_CGMFtk():
    try:
        import CGMFtk
        print("CGMFtk already installed")
    except ImportError as e:
        print("CGMFtk has not been installed yet, installing now...")
            # Change directory to tools
        os.chdir(os.path.join(os.path.dirname(__file__), 'tools'))

        # Use subprocess to run a shell command to build and install CGMFtk
        subprocess.check_call(['python', 'setup.py', 'install'])

        os.chdir(os.path.join(os.path.dirname(__file__)))

# first install CGMFtk
#install_CGMFtk()

# then run setup
setup(
    name="pyCGMF",
    version="1.0",
    description="Python drivers to run histories of CGMF in parallel using MPI and analyze results in Python",
    author="Kyle Beyer",
    author_email="beykyle@umich.edu",
    zip_safe=False,
    ext_modules=[CMakeExtension("pyCGMF")],
    cmdclass={"build_ext": CMakeBuild},
    extras_require={"test": ["pytest>=6.0"]},
    python_requires=">=3.7",
)

