from Cython.Build import cythonize
from distutils.core import setup, Extension
#from setuptools import setup, Extension
import Cython
import multiprocessing
import os
import shutil

print("Using Cython %s" % Cython.__version__)

def configure_google_hashmap():
    # TODO: In temp dir or something
    print("Configuring Google hash map")
    if os.system("3rd-party/sparsehash/configure") != 0:
        raise RuntimeError("Error configuring Google hash map")
    shutil.copy("src/config.h", "cpp/sparsehash/internal/sparseconfig.h")

def slurp(filename):
    with open(filename, "rt") as f:
        return f.read()

configure_google_hashmap()

extensions = [
    Extension("_arv", [
            "cpp/dnatraits.cpp",
            "cpp/file.cpp",
            "cpp/fileptr.cpp",
            "cpp/filesize.cpp",
            "cpp/mmap.cpp",
            "cpp/parse_file.cpp",
            "cython/_arv.pyx",
        ],
        language="c++",
        include_dirs=["cpp"],
        extra_compile_args=["-std=c++11", "-g0"], # gcc specific
    ),
]

keywords = [
    "23andMe",
    "bio",
    "biology",
    "biopython",
    "disease",
    "DNA",
    "gene",
    "genome",
    "health",
    "protein",
    "RNA",
    "RSID",
    "SNP",
]

setup(
    name="arv",
    packages=["arv"],
    package_dir={"arv": "python/arv"},
    version="0.1",
    description="A fast 23andMe raw genome file parser",
    author="Christian Stigen Larsen",
    author_email="csl@csl.name",
    url="https://github.com/cslarsen/arv",
    license="https://www.gnu.org/licenses/gpl-3.0.html",
    long_description=slurp("README.rst"),
    keywords=keywords,
    platforms=["unix", "linux", "osx"],
    ext_modules=cythonize(extensions, nthreads=multiprocessing.cpu_count()),
)
