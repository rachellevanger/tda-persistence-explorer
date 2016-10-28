import IPython
import json
import random
import string


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
