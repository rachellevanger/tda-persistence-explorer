# TDA Persistence Explorer

A data exploration tool for Topological Data Analysis (TDA). This repository contains a D3-based Python module and a sample Jupyter notebook suitable for analyzing persistence diagrams generated from 2D images, and is especially useful for studying time series of images (i.e. video). 


Contributors:
Jacek Cyranka,
Shaun Harker,
Rachel Levanger


# Docker Installation

We recommend trying our Docker installation of the app, unless you are a Python ninja, in which case feel free to follow the instructions under the Local Installation heading, below.

To run tda-persistence-explorer via Docker, you'll first need to [download and install Docker](https://www.docker.com/) if you don't already have it. Once you have installed Docker, run the following command in a terminal to pull down our docker image:

`docker pull rachellevanger/tda-persistence-explorer`

**For MAC users:**

Make sure you have the `/Users` folder shared. From the Docker menu, choose Preferences... and then go to the File Sharing tab. This folder is shared by default, so it should already be in the list. If it isn't, add it and then restart Docker.

*Try out the Tutorial first:*

To try out tda-persistence-explorer for the first time, we recommend running the pre-packaged tutorial. To start it up, copy/paste the following command into a terminal:

`docker run -it -p 8888:8888 rachellevanger/tda-persistence-explorer sh -c "cd tda-persistence-explorer; jupyter notebook --ip=0.0.0.0 --no-browser"`

The output from this command will include something that looks like this:

```
Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://0.0.0.0:8888/?token=2b87eb399551840a2d418489d93e68664e4dfb2859193379
```
Copy/paste the link, including the token, and this will open up a Jupyter notebook interface. Navigate to `/docs/Tutorial.ipynb` and then follow the instructions provided in the notebook.

*To run your own notebooks on your local filesystem:*

Once you are familiar with the Tutorial, you're ready to start working with your own data! You will need to tell the Docker container to mount your local `/Users` directory in order to browse to files you might have in your `Documents` folder. (You'll need to specify a different root direcotry if you wish to work with data outside of the `/Users` directory, and make sure it is shared with Docker.)

We suggest running the following command to start getting familiar with running tda-persistence-explorer through this Docker container:

`docker run -it -p 8888:8888 -v /Users:/Users rachellevanger/tda-persistence-explorer sh -c "jupyter notebook --ip=0.0.0.0 --no-browser"`

Copy/paste the link with the provided token as you did for the Tutorial. The `/Users` directory should be available from the list of directories. Browse to your `.ipynb` file and click on it to run the notebook. It shoudl now have access to the PersistenceExplorer installation in the Docker container.



# Local Installation

Type the following from the root of repository in order to install the program:

    ./install.sh

The program has the following dependencies:

* C++11 compiler (any modern gcc or clang will do)
* CImg
* PHAT

## Troubleshooting

The installer attempts to download the CImg and PHAT dependencies. If `wget` or `git` are not available on your system, this step will fail and the installation will abort. `wget` and `git` are important utilities for downloading files from the internet from the command line.

On Mac OS X, `wget` is best obtained with Homebrew.

## Try it out

Once you've installed the app, type the following in the command line to start up a jupyter notebook:

     jupyter notebook

(If you don't have jupyter notebooks installed, go here: <http://jupyter.readthedocs.io/en/latest/install.html>.)

Once your local notebook server is running, navigate to `/doc/Tutorial.ipynb` and start up the notebook. The Tutorial will walk you through a sample dataset provided with this installation.

# Maintainer Notes

The `ImagePersistence` program is built and placed in `/source/PersistenceExplorer/bin`. The Python package installer then copies this executable onto path so it can be used.

Python packaging: <http://python-packaging.readthedocs.io/en/latest/everything.html>



