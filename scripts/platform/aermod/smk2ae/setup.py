from setuptools import setup, find_packages

setup(
    name="smk2ae",
    version="0.2",
    packages=['smk2ae'],
    package_dir = {'smk2ae': 'smk2ae'},
    python_requires='>3.5',
    scripts=['scripts/cmv_smk2ae.py','scripts/aermod.py'],
    setup_requires=['pandas>=0.20','numpy>=1.12','pyproj>=1.9','GDAL>=2.2'],
    author_email='beidler.james@epa.gov'
)
