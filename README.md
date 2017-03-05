# TDA Persistence Explorer

A data exploration tool for Topological Data Analysis (TDA). This repository contains a D3-based Python module and a sample Jupyter notebook suitable for analyzing persistence diagrams generated from 2D images, and is especially useful for studying time series of images (i.e. video). 


Contributors:
Jacek Cyranka,
Shaun Harker,
Rachel Levanger

We provide two sets of installation instructions:

* [Docker Installation](#docker-installation) - Simple
* [Local Installation](#local-installation) - Advanced


# Docker Installation

For your convenience, we provide a Docker image that comes with everything in this repository (and all of the dependencies) pre-installed. 

**System Requirements:**

* [Docker](https://www.docker.com/)
* Approx 1.5GB of disk space (for the image)
* 4GB of RAM (required by Docker)

To run tda-persistence-explorer via Docker, you'll first need to download and install Docker if you don't already have it. Once you have installed Docker and have verified the Docker service is running, run the following command in a terminal to pull down our docker image:

`docker pull rachellevanger/tda-persistence-explorer`

After the image is pulled, you can verify it is avaialble by running `docker images`.

## Running the Tutorial Jupyter notebook via a Docker container

To run Jupyter notebook from the Docker container, a port will need to be exposed on your computer. To secure your connection to the notebook server, we suggest using the following command to start up the notebook server.

```
docker run -d -e GEN_CERT=yes -p 8888:8888 rachellevanger/tda-persistence-explorer /bin/bash -c "cd tda-persistence-explorer; start-notebook.sh" && sleep 2 && docker logs $(docker ps -l -q) 2>&1 | grep https://localhost:
```

Copy the URL `https://localhost:8888/[very long token id]` from the output and paste it into your browser. Your browser will warn you that this site is not trusted. The option `-e GEN_CERT=yes` in the above command instructed the container to generate a self-signed SSL certificate and configured the Jupyter notebook server to accept HTTPS connections. Typically, sites with self-signed certificates should not be trusted, but since you are the one who created this site, this is okay. Ignore the scary messages and continue.

To try out tda-persistence-explorer for the first time, browse to `doc/Tutorial.ipynb`. Follow the instructions provided in the notebook to explore a sample dataset.

## Running your own notebooks on your local filesystem via a Docker container

Once you are familiar with the Tutorial, you're ready to start working with your own data! You will need to tell the Docker container to mount to your local data directory in order to access data on your local computer. Note that this gives anyone with your public IP address and the Jupyter notebook token access to your data directory, which is why it's best practice to use HTTPS.

Docker should already have access to common system folders (e.g. `/Users` for Macs), since this is the default setup. To check that Docker has access to the data directory you will need to access, from the Docker menu choose Preferences... and then go to the File Sharing tab. If your directory isn't contained in a directory in the list, add it and then restart Docker.

To run PersistenceExplorer from a Jupyter notebook located outside of the Docker container, copy/paste the following into a terminal, modifying the path `/path/to/my/local/work/directory` so that it points to the directory with your working data in it:

```
docker run -d -e GEN_CERT=yes -v /path/to/my/local/work/directory:/home/jovyan/work -p 8888:8888 rachellevanger/tda-persistence-explorer start-notebook.sh && sleep 2 && docker logs $(docker ps -l -q) 2>&1 | grep https://localhost:
```

Copy the URL `https://localhost:8888/[very long token id]` from the output and paste it into your browser. Since this command generates a self-signed SSL certificate, your browser will warn you that the link is untrusted. Ignore the browser warnings in order to run the notebook. The default directory from the Jupyter noteook will be the local folder that you mounted. Browse to your `.ipynb` file and click on it to run the notebook. It should now be running via the installation of the tda-persistence-explorer app in the Docker container.

Note, you can also change the port number (here we choose 8888) in the event you want multiple containers running for exploring different datasets simultaneously (provided your computer has enough resources to do so).


## Keeping your Docker environment clean

Each time you run `docker run` on the commandline, a new Docker container is instantiated. Learn how to keep your docker environment clean by reading [this helpful guide](https://www.digitalocean.com/community/tutorials/how-to-remove-docker-images-containers-and-volumes).

It is easiest (and most secure) to stop and remove the Docker container when you are done using it by running the following commands, which will stop and then remove all docker containers.

```
docker stop $(docker ps -a -q)
docker rm $(docker ps -a -q)
```


# Local Installation

After cloning the repository to a local direcotry, type the following from the root of repository in order to install the program:

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

The Docker image is based on the [jupyter/minimal-notebook](https://github.com/jupyter/docker-stacks/tree/master/minimal-notebook) image.

