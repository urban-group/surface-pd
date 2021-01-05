from setuptools import setup, find_packages
import os
import glob

__author__ = "Xinhao Li"
__email__ = "xl2778@columbia.edu"

here = os.path.abspath(os.path.dirname(__file__))
package_name = 'surface_pd'
package_description = 'Generation of surface phase diagrams'

# Get the long description from the README file
with open(os.path.join(here, 'README.md')) as fp:
    long_description = fp.read()

# Get version number from the VERSION file
with open(os.path.join(here, package_name, 'VERSION')) as fp:
    version = fp.read().strip()

setup(
    name=package_name,
    version=version,
    description=package_description,
    long_description=long_description,
    # url='',
    author=__author__,
    author_email=__email__,
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    keywords=['materials science', 'phase diagrams'],
    packages=find_packages(exclude=['tests']),
    install_requires=['numpy', 'pymatgen'],
    # package_data={
    #     'sample': [''],
    # },
    scripts=glob.glob(os.path.join("scripts", "*.py"))
)
