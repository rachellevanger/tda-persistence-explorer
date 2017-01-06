# TDA D3 Explorer

A data exploration tool for Topological Data Analysis (TDA). This repository contains a D3-based Python module and a sample Jupyter notebook suitable for analyzing persistence diagrams generated from 2D images, and is especially useful for studying time series of images (i.e. video). 


Contributors:
Jacek Cyranka,
Shaun Harker,
Rachel Levanger


# Installation

Type the following from the root of repository in order to install the program:

    ./install.sh

The program has the following dependencies:

* C++11 compiler (any modern gcc or clang will do)
* CImg
* PHAT

## Troubleshooting

The installer attempts to download the CImg and PHAT dependencies. If `wget` or `git` are not available on your system, this step will fail and the installation will abort. `wget` and `git` are important utilities for downloading files from the internet from the command line.

On Mac OS X, `wget` is best obtained with Homebrew.

# Try it out

Once you've installed the app, type the following in the command line to start up a jupyter notebook:

     jupyter notebook

(If you don't have jupyter notebooks installed, go here: <http://jupyter.readthedocs.io/en/latest/install.html>.)

Once your local notebook server is running, navigate to `/doc/Tutorial.ipynb` and start up the notebook. The Tutorial will walk you through a sample dataset provided with this installation.

# Maintainer Notes

The `ImagePersistence` program is built and placed in `/source/PersistenceExplorer/bin`. The Python package installer then copies this executable onto path so it can be used.

Python packaging: <http://python-packaging.readthedocs.io/en/latest/everything.html>



