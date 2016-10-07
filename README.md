# TDA D3 Explorer

A D3 exploration tool for Topological Data Analysis (TDA). This repository contains a D3 .html page suitable for analyzing persistence diagrams generated from 2D images, and is especially useful for studying time series of images (i.e. video). 


Contributors:
Jacek Cyranka,
Shaun Harker,
Rachel Levanger


To try it out, after cloning the repository, run the following from the command line in the repository directory (assuming you have Python installed):

`python -m SimpleHTTPServer 8000`

If this doesn't work, use your favorite local HTTP server. Then navigate to the following url:

[http://localhost:8000/explorer.html?imgStart=1&trail=5&listGens=0&selector=1](http://localhost:8000/explorer.html?imgStart=1&trail=5&listGens=0&selector=1)

The arguments in the url give a few parameters used to control the web tool:
* imgStart = This is the starting index of the tool (1 to 20 for the supplied test data)
* trail = This is the number of subsequent frames to overlay on the persistence diagrams (e.g. 0 to show just a single image at a time, 5 to show the current frame plus five more)
* listGens = 0 to suppress a generator listing, 1 to show it
* selector = Controls which persistence diagram is active for the selection tool. Choose 1 for top left, 2 for top right, 3 for bottom left, and 4 for bottom right.

Once the page loads, on the left you'll see a beautiful picture of simulated convection flows (simulations generously provided by the [Paul Research Group](http://www.me.vt.edu/mpaul/) at Virginia Tech). On the right you'll see both sublevel and superlevel persistence diagrams. Click and drag to circle points on the selected persistence diagram, and when you release, the images will play forward in time with generators corresponding to the encircled persistence points overlaid in cyan. The open circle marks the location of the birth critical cell and the closed circle marks the location of the death critical cell.

To use the tool with your own data, first generate the persistence diagrams with the code contained in `phat_persistence_from_image/`, following the instructions in the README on how to compile and run the software. After your data is generated, you'll have to read the `explorer.html` file and change a whole lotta stuff in order to make it work. We're working on this to make it easeir to use, we promise!!

