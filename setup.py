import setuptools
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
        name="colltree",
        version="0.0.1",
        author="Jennifer Scora",
        author_email="jennifer.scora@mail.utoronto.ca",
        description="A package to find the collision history and interactively plot it",
        long_description=long_description,
        long_description_content_type="text/markdown",
        classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python'
        ],
        url='https://github.com/jscora/colltree',
        license='MIT',
        packages=["colltree"],
        package_data={'colltree':['assets/dash_components.css']},
        install_requires=["numpy","astropy","pandas","jupyter_dash","plotly","tqdm"]
        )

