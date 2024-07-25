from setuptools import setup, find_packages

setup(
    name='meshpyProcessing',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'requests',
        'owslib'
    ],
    description='Python package for extracting and processing hydrometric data.',
    author='Fuad Yassin',
    author_email='fuad.yassin@usask.ca',
    url='https://github.com/fuadyassin/meshpyProcessing',
)
