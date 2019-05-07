import os

from setuptools import setup


def get_requirements():
    """Get requirements from requirements.txt"""
    with open("requirements.txt") as handle:
        return [x.strip() for x in handle]


def get_version():
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(os.path.join("peddy", "__init__.py"), "r") as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
             if isinstance(node, ast.Assign)
             and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")

setup(
    version=get_version(),
    name='peddy',
    description="pleasingly pythonic pedigree manipulation",
    packages=['peddy', 'peddy.tests'],
    long_description=open('README.md').read(),
    author="Brent Pedersen",
    author_email="bpederse@gmail.com",
    install_requires=get_requirements(),
    zip_safe=False,
    test_suite='nose.collector',
    # To provide executable scripts, use entry points in preference to the
       # "scripts" keyword. Entry points provide cross-platform support and
       # allow pip to create the appropriate form of executable for the
       # target platform.
    entry_points=dict(
        console_scripts=[
            'peddy = peddy.__main__:cli',
        ],
    ),
    include_package_data=True,
    tests_require='nose',
    classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
)
