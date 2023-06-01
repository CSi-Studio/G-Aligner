import os
import sys
import glob
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import setuptools

__version__ = '0.9.4'
requirements_list = ['pybind11>=2.4', 'numpy']

maindir = os.path.join(".", "obiwarp")
libdir = os.path.join(maindir, "lib")
source_files = ['py_obiwarp.cc']
file_types = ['*.c', '*.cc', '*.cpp']
exclude_patterns = ['TEMPLATE', '2']
libraries = []
extra_objects = []

for file_type in file_types:
    for f in glob.glob(os.path.join(libdir, file_type)):
        exclude = False
        for pattern in exclude_patterns:
            if pattern in f:
                exclude = True
                break
        if not exclude:
            source_files += [f]

ext_modules = [
    Extension(
        'py_obiwarp',
        source_files,
    ),
]

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    #The c++14 is prefered over c++11 (when it is available).
    # This somehow can fail on a Mac with clang
    #"""
    if has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc', '/openmp', '/O2'],
        'unix': ['-O3', '-march=native', '-std=c99'],
        #'unix': ['-O0', '-march=native', '-g'],
    }
    link_opts = {
        'unix': [],
        'msvc': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        link_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    else:
        c_opts['unix'].append("-fopenmp")
        link_opts['unix'].extend(['-fopenmp', '-pthread'])

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())

        # extend include dirs here (don't assume numpy/pybind11 are installed when first run, since
        # pip could have installed them as part of executing this script
        import pybind11
        import numpy as np
        for ext in self.extensions:
            ext.extra_compile_args.extend(opts)
            ext.extra_link_args.extend(self.link_opts.get(ct, []))
            ext.include_dirs.extend([
                # Path to pybind11 headers
                pybind11.get_include(),
                pybind11.get_include(True),

                # Path to numpy headers
                np.get_include()
            ])

        build_ext.build_extensions(self)


setup(
    name='py_obiwarp',
    version='0.9.4',
    author="CSi-Studio",
    author_email="csi@csibio.net",
    maintainer="Ruimin Wang",
    license="The MIT License",
    description='Python bindings for OBI-Warp',
    long_description='Python bindings for OBI-Warp algorithm, which was also used in XCMS for mass spectrometry data alignment.',
    keywords="ObiWarp, Alignment, Warping",
    packages=find_packages(),
    ext_modules=ext_modules,
    install_requires=requirements_list,
    setup_requires=requirements_list,
    cmdclass={'build_ext': BuildExt},
    test_suite="tests",
    zip_safe=False,
    python_requires=">=3.6"
)
