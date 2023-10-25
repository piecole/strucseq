from setuptools import setup, find_packages

setup(
    name='strucsec',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'requests',
        'beautifulsoup4',
        'pandas',
        'lxml',
        'numpy'
    ],
    author='Pierre Coleman',
    description='A package for fetching and processing data between protein sequences and structure',
    url='https://github.com/piecole/strucseq'
)
