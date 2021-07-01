import setuptools

setuptools.setup(
    name="pycatkin", 
    version="0.1.0",
    author="Astrid Boje",
    description="Python Catalysis Kinetics Toolset",
    packages=["pycatkin"],
    python_requires=">=3.7",
    install_requires=["numpy",
                      "scipy",
                      "matplotlib",
                      "ase",
                      "pandas",
                      "sphinx",
                      "sphinx_rtd_theme"
                      ]
)
