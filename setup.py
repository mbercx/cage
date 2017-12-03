from setuptools import setup, find_packages

setup(
    name='cage',
    version='0.1',
    packages=find_packages(exclude=['docs']),
    install_requires=[
        'matplotlib',
        'scipy',
        'numpy',
        'click',
        'pymatgen'
    ],
    entry_points='''
        [console_scripts]
        cage=cage.cli.cli:main
    ''',
)