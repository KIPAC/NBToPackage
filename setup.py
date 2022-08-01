from setuptools import setup

setup(
    name="piff_syst_plots",
    author="Eric Charles",
    author_email="echarles@slac.stanford.edu",
    url="https://github.com/KIPAC/NBToPackage",
    packages=["piff_syst_plots"],
    description="Plots for PIFF systematics studies",
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
    long_description=open("README.md").read(),
    package_data={"": ["README.md", "LICENSE", "*.npy"]},
    use_scm_version={"write_to": "piff_syst_plots/_version.py"},
    include_package_data=True,
    scripts=["bin/piff_syst_plots"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    install_requires=[
        "matplotlib",
        "numpy",
    ],
)
