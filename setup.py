from setuptools import setup, find_packages

setup(
    name = "cirrus-ngs",
    version = "1.0.0a1",
    description = "Cloud-optimized primary analysis pipelines for RNA-seq, miRNA-seq, ChIP-seq, and variant calling in whole-genome/whole-exome DNA-seq",
    long_description_content_type = "text/markdown",
    long_description = open("README.md").read(),
    url = "https://github.com/ucsd-ccbb/cirrus-ngs",
    license = "MIT",
    author = "Mustafa Guler",
    author_email = "mguler@ucsd.edu",
    classifiers = 
    [
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    keywords = "bioinformatics, NGS, jupyter",
    install_requires = 
    [
        "paramiko",
        "pyyaml",
        "jupyter",
        "notebook",
        "scp",
        "cfncluster"
    ],
    python_requires = ">=3.5",
    packages = find_packages(exclude = ("Pipelines*", "server*", "tests*")),
    maintainer = "Mustafa Guler",
    maintainer_email = "mguler@ucsd.edu",
    entry_points = 
    {
        "console_scripts": ["cirrus-ngs = cirrusngs.cli:main"]
    },
    include_package_data=True
)
