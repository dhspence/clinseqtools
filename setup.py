from setuptools import setup, find_packages

setup(
    name='clinseqtools',  # Name of your package
    version='0.1',        # Version of your module
    packages=find_packages(),
    install_requires=[
        'pandas','pysam'
    ],
)
