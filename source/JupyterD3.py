import IPython
import json
import random
import string
import os
import subprocess

def Execute(command_string, html_string):
  """
  Execute(command_string, html_string)
    Submits the command to be executed as javascript along with html_string
    D3 and PersistenceExplorer are preloaded if they aren't already
  """
  with open("PersistenceExplorer.css") as f:
    stylesheet = "<style>" + f.read() + "</style>"
  output = stylesheet + """
    <script>
    var command = function() { 
    """ + command_string + """ };
    function LoadSource(src, tailcall) {
      var elements = document.querySelectorAll("script[src='"+src+"']");
      if ( elements.length == 0 ) {
        var element = document.createElement("script");
        element.src = src;
        document.body.appendChild(element);
        element.onload = tailcall;
      } else {
        tailcall ();
      }
    };
    LoadSource("//d3js.org/d3.v3.min.js", function() {
      LoadSource( "./PersistenceExplorer.js", function() {  
          command();
      })
    });
    </script>
    """ + html_string
  return IPython.display.HTML(output)


def PersistenceExplorer( imagefiles, persistencefiles, frames, dimension ):
  # Create a unique prefix identifier to prevent confused ids 
  prefix = ''.join(random.SystemRandom().choice(string.ascii_uppercase) for _ in range(16))
  # Substitute passed parameters into call
  command_string = "PersistenceExplorer(_imagefiles, _persistencefiles, _frames, _dimension, _divs);"
  command_string = command_string.replace('_imagefiles', json.dumps(imagefiles))
  command_string = command_string.replace('_persistencefiles', json.dumps(persistencefiles))
  command_string = command_string.replace('_frames', json.dumps(frames))
  command_string = command_string.replace('_dimension', json.dumps(dimension))
  command_string = command_string.replace('_divs', 
    json.dumps({'divImg' : '#' + prefix + 'divImg', 
                'divPD' : '#' + prefix + 'divPD',
                'divSlide' : '#' + prefix + 'divSlide'}))
  # Create HTML string to write app into
  html_string = """
<div id="divImg" class="PersistenceExplorerPanel"></div>
<div id="divPD" class="PersistenceExplorerPanel"></div>
<br />
<p></p>
<div id="divSlide" class="PersistenceExplorerPanel" ></div>
<p></p>
"""
  html_string = html_string.replace("divImg", prefix + "divImg" )
  html_string = html_string.replace("divPD", prefix + "divPD" )
  html_string = html_string.replace("divSlide", prefix + "divSlide" )

  # Display the HTML
  return Execute(command_string, html_string)

def ProcessImageListWithPHAT( list_of_image_filenames, list_of_output_filenames, filtration_type ):
  """
  Iterate through images, compute persistence results, and store results.
    list_of_image_filenames: a list of image files
    list_of_output_filenames: a list of files to save corresponding persistence results in
    filtration_type: either "sub" or "super" to indicate to obtain persistence results for either
                     sublevel or superlevel set filtrations
  """
  # Run commands in parallel
  processes = [subprocess.Popen(["./2d_pic_to_bd_matrix", infile, outfile, filtration_type]) for infile, outfile in zip(list_of_image_filenames, list_of_output_filenames) ]
  # Block until processing complete
  exitcodes = [p.wait() for p in processes]


def ProcessImageFolderWithPHAT( image_foldername, sub_foldername=None, sup_foldername=None ):
  """
  Compute sublevel and superlevel persistence diagrams for all images in a folder.
    image_foldername: path to images
    sub_foldername: path to put sublevel persistence results (defaults to image_foldername/sub)
    sub_foldername: path to put sublevel persistence results (defaults to image_foldername/sup)
  """
  list_of_image_filenames = [filename for filename in os.listdir(image_foldername) if os.path.isfile(os.path.join(image_foldername,filename)) and not filename.endswith('.md') and not filename.endswith('.txt') ]
  list_of_output_filenames = [os.path.splitext(filename)[0] for filename in list_of_image_filenames]
  if sub_foldername is None:
    sub_foldername = image_foldername + "/pd_sub"
  if sup_foldername is None:
    sup_foldername = image_foldername + "/pd_sup"
  if not os.path.exists(sub_foldername):
    os.makedirs(sub_foldername)
  if not os.path.exists(sup_foldername):
    os.makedirs(sup_foldername)
  list_of_image_filenames = [ os.path.join(image_foldername,filename) for filename in list_of_image_filenames ]
  list_of_sub_output_filenames = [ os.path.join(sub_foldername,filename) for filename in list_of_output_filenames ]
  list_of_sup_output_filenames = [ os.path.join(sup_foldername,filename) for filename in list_of_output_filenames ]
  ProcessImageListWithPHAT(list_of_image_filenames, list_of_sub_output_filenames, 'sub')
  ProcessImageListWithPHAT(list_of_image_filenames, list_of_sup_output_filenames, 'super')
