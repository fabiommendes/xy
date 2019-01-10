import os
from pathlib import Path

from setuptools import setup, find_packages
from setuptools.extension import Extension

BASE = Path(__file__).parent
INCLUDES = [str(BASE / 'include')]

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext, Extension
except ImportError:
    setup_kwds = {}
    cythonize = (lambda x, **kwargs: x)
    print('Warning: Cython not found!')
else:
    setup_kwds = {'cmdclass': {'build_ext': build_ext}}


def extension(name: str, compile_args=(), link_args=(), include_dirs=(),
              libraries=(), **kwargs):
    """
    Standard Cython extension.
    """
    path = os.path.sep.join(['src', *name.split('.')]) + '.pyx'
    include_dirs = ['include', *include_dirs]
    return Extension(name,
                     [path],
                     extra_compile_args=compile_args,
                     extra_link_args=link_args,
                     include_dirs=include_dirs,
                     libraries=libraries,
                     **kwargs)


def get_extensions():
    return cythonize([
        extension('xy.linalg2d.vector_2d', language='c++'),
    ], include_path=['include'])


# Meta information
version = '0.1.0b0'

setup(
    # Basic info
    name='xy-math',
    version=version,
    author='Fábio Macêdo Mendes',
    author_email='fabiomacedomendes@gmail.com',
    url='http://github.com/fabiommendes/smallvectors/',
    description='Linear algebra objects and shapes in small dimensions.',
    long_description=open('README.rst').read(),

    # Classifiers (see https://pypi.python.org/pypi?%3Aaction=list_classifiers)
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Software Development :: Libraries',
    ],

    # Packages and dependencies
    package_dir={'': 'src'},
    packages=find_packages('src'),
    install_requires=[
    ],
    extras_require={
        'dev': [
            'pytest',
            'manuel',
            'sphinx',
            'invoke',
        ],
    },
    package_data={'euclid.includes': ['*.pxi'],
                  'euclid.linalg2d': ['*.pxd', '*.pyx'],
                  'euclid.shapes2d': ['*.pxd', '*.pyx']},

    # Other configurations
    platforms='any',
    zip_safe=False,
    ext_modules=get_extensions(),
    **setup_kwds,
)
