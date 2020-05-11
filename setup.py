import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()


setuptools.setup(
    name = "adpred",
    version="0.0.2",
    author="Ariel Erijman",
    author_email="aerijman@fredhutch.org, aerijman@neb.com",
    description="Prediction of Transcription Activation Domains from protein sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FredHutch/adpred-pkg",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
