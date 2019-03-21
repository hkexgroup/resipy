# to build run : python3 setup.py sdist bdist_wheel
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyR2",
    version="1.1.1",
    author="HKEx",
    description="API for ERT inversion (DC/IP)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hkex/pyr2",
    packages=setuptools.find_packages(),
    install_requires=['numpy','matplotlib','pandas','scipy','statsmodels'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: OS Independent",
    ],
)
