from setuptools import setup, find_packages
from pathlib import Path

VERSION = "1.0.2"
DESCRIPTION = "Generalized Explicit Runge-Kutta"

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="gerk",
    version=VERSION,
    author="Yasser Naji",
    author_email="yfnaji@gmail.com",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=["numpy"],
    keywords=[
        "runge-kutta", 
        "runge", 
        "kutta", 
        "numerical", 
        "integration", 
        "approximation"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)