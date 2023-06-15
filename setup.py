from setuptools import setup, find_packages
import codecs
import os

VERSION = "0.0.1"
DESCRIPTION = "Generalized Explicit Runge-Kutta"

setup(
    name="gerk",
    version=VERSION,
    author="Yasser Naji",
    author_email="yfnaji@gmail.com",
    description=DESCRIPTION,
    packages=find_packages(),
    install_requires=["matplotlib"],
    keywords=["runge-kutta", "runge", "kutta", "numerical", "integration", "approximation"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)