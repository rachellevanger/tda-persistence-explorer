import IPython
import json

def LoadD3():
	return IPython.display.HTML("""
		<script src="//d3js.org/d3.v3.min.js"></script>
		D3 Loaded.
	""")

def PersistenceExplorer( imagefiles, persistencefiles, frames, dimension ):
  # Substitute passed parameters into call
  command_string = "<script> PersistenceExplorer(_imagefiles, _persistencefiles, _frames, _dimension); </script>"
  command_string = command_string.replace('_imagefiles', json.dumps(imagefiles));
  command_string = command_string.replace('_persistencefiles', json.dumps(persistencefiles));
  command_string = command_string.replace('_frames', json.dumps(frames));
  command_string = command_string.replace('_dimension', json.dumps(dimension));
  # Load style sheet and javascript
  with open("PersistenceExplorer.css") as f:
    stylesheet = "<style>" + f.read() + "</style>"
  with open("PersistenceExplorer.js") as f:
    javascript = "<script>" + f.read() + "</script>"
  # Div setup
  html_string = """
<div id="imgContainer"></div>
<div id="divImg" class="PersistenceExplorerPanel"></div>
<div id="divPD" class="PersistenceExplorerPanel"></div>
<br />
<p></p>
<div id="slide" class="PersistenceExplorerPanel" ></div>
<p></p>
"""
  # Display the HTML
  return IPython.display.HTML(stylesheet + html_string + javascript + command_string);
