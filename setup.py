from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
from glob import glob

__version__ = "0.0.1"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)


ext_modules = [
    Pybind11Extension(
        "control_and_sim_algorithms",
        #sorted(glob("src/*.cpp")),
        ["src/main.cpp", "src/solvers.cpp", "src/systems.cpp"], # Important: include all cpp files
        # Example: passing in the version to the compiled code
        define_macros=[("VERSION_INFO", __version__)],
        cxx_std=17,
    ),
]

setup(
    name="control_and_sim_algorithms",
    version=__version__,
    author="Roman Frels",
    description="A collection of control and simulation algorithms, implemented in C++ and Python via pybind11",
    long_description="",
    package_dir={"": "."},
    packages=[],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.11",
)
