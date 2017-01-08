// PersistenceExplorer.js
// MIT LICENSE. Rachel Levanger 2016.
// Contributors: Rachel Levanger (author) 
//               Shaun Harker (contributor)

// Requires: d3



// String.prototype.format
//   Source: http://stackoverflow.com/a/4256130/1467617
//   Effects: Acts like sprintf for Javascript.
//   Example: 'The {0} is dead. Don\'t code {0}. Code {1} that is open source!'.format('ASP', 'JS');
String.prototype.format = function() {
    var formatted = this;
    for (var i = 0; i < arguments.length; i++) {
        var regexp = new RegExp('\\{'+i+'\\}', 'gi');
        formatted = formatted.replace(regexp, arguments[i]);
    }
    return formatted;
};

// loadImages(_imagefilenames, _frames, _divs)
//   Inputs:
//     _imagefilenames : a list of filenames pointing to image data
//     _frames : a list of frames (corresponding to _imagefilenames indexing) of interest
//     _divs : the CSS selectors indicating the divs to draw app in
//   Effects:
//     Loads the images _imagefilenames[_frames[i]] for 0 <= i < _frames.length.
//     More specifically, it creates img tags with src files in _imagefilenames 
//     for the indicated _frames.
//     It expects there is a "div" named "divImg"
//     It creates an <img> tag inside this div with attribute "frame_number"  
//     for each frame_number in _frames
//     and sets its "src" attribute to the filename _imagefilenames[frame_number] 
function loadImages(_imagefilenames, _frames, _divs) {
  _frames.forEach(function(frame_number){
      d3.select(_divs.divImg)
          .append('img')
            .attr('src', _imagefilenames[frame_number])
            .attr('frame_number', frame_number)
            .attr('class', 'img')
            .attr("style","display: none;");
  });
}

// loadPersistenceData(_files, _frames)
//   Inputs:
//     _files : a list of filenames pointing to persistence data.
//              Each file is a .csv file which contains fields for 
//              "birth" "death" "b_x" "b_y" "d_x" "d_y" and "dim"
//     _frames : a list of frames (corresponding to _files indexing) of interest
//     _scale : a scaling factor for the birth and death critical cell locations
//   Outputs:
//     Loads the files _files[_frames[i]] for 0 <= i < _frames.length, and returns
//     a promise-wrapped list of objects with "birth" "death" "time" and "selected" fields.
//     The "time" field gives the index of the file in which the data was loaded.
//     The "selected" field is initialized as "false"
function loadPersistenceData(_files, _frames, _scale) {
  function PromiseToLoadSingleFile(frame_number){
    return new Promise(function(resolve, reject){
        var img = new Image()
        d3.csv(_files[frame_number], function(error, dataset) {
          dataset.forEach(function(d) {
            d.birth = parseInt(d.birth);
            d.death = parseInt(d.death);
            d.b_x = parseInt(parseFloat(d.b_x)*_scale);
            d.b_y = parseInt(parseFloat(d.b_y)*_scale);
            d.d_x = parseInt(parseFloat(d.d_x)*_scale);
            d.d_y = parseInt(parseFloat(d.d_y)*_scale);
            d.time = frame_number;
            d.selected = false;
          });
          resolve(dataset);
        });
      });
  };
  var promise_to_load_all_files = Promise.all(_frames.map(PromiseToLoadSingleFile));
  var concatenate_lists_together = function(lists) {return [].concat.apply([], lists);};
  return promise_to_load_all_files.then(concatenate_lists_together);
}

// isFeatureSelected(_d,_extent)
//   Inputs:
//     _d : persistence point object with fields b_x, b_y, d_x, and d_y
//     _extent : a rectangle stored as [ [ xmin, ymin ], [xmax, ymax] ]
//   Output:
//     Returns true if either [b_x, b_y] or [d_x, d_y] is in the rectangle _extent
function isFeatureSelected(_d, _extent) {
  return ( ( _extent[0][0] <= _d.b_x ) && ( _d.b_x <= (_extent[1][0]) ) && ( _extent[0][1] <= _d.b_y ) && ( _d.b_y <= _extent[1][1] )
    || ( _extent[0][0] <= _d.d_x ) && ( _d.d_x <= (_extent[1][0]) ) && ( _extent[0][1] <= _d.d_y ) && ( _d.d_y <= _extent[1][1] ) );
}

// annotateImageWithFeatures(_data, _frame_number, _svg, _dimension )
//   Input: 
//     _data : list of persistence points to be plotted
//     _frame_number : frame number of interest ( frame number, not frame index i.e. _frame[frame_index] == _frame_number )
//     _svg : SVG containing image annotations
//     _dimension : an integer indicating dimension of persistence diagram of interest
//   Effects:
//     Replaces annotations in _svg with the ones corresponding to 
//     the selected ones in _data for frame _frame_number. The annotations
//     are line segments with open circles at the birth feature endpoints
//     and closed circles at the death feature endpoints.
function annotateImageWithFeatures( _data, _frame_number, _svg, _dimension ){

  // Remove previous drawings done in previous calls to "imageUpdate"
  _svg.selectAll(".dot_img")
    .remove();

  // Create unfilled circles of radius 3 for all selected birth features
  // in frame _frame_number with given dimension
  _svg.selectAll(".dot")
      .data(_data.data().filter(function(d) {return d.time == _frame_number & d.dim == _dimension & d.selected===true}))
      .enter()
    .append("circle")
      .attr("class", "dot_img")
      .attr("r", 3)
      .attr("cx", function(d) { return d.b_x; })
      .attr("cy", function(d) { return d.b_y; })
      .style("stroke","cyan")
      .style("fill", 'none');

  // Create filled circles of radius 3 for all selected death features
  // in frame _frame_number with given dimension
  _svg.selectAll(".dot")
      .data(_data.data().filter(function(d) {return d.time == _frame_number & d.dim == _dimension & d.selected===true}))
      .enter()
    .append("circle")
      .attr("class", "dot_img")
      .attr("r", 3)
      .attr("cx", function(d) { return d.d_x; })
      .attr("cy", function(d) { return d.d_y; })
      .style("stroke","cyan")
      .style("fill", "cyan");

  // Create line segments connecting birth and death features
  // in frame _frame_number with given dimension
  _svg.selectAll(".dot")
      .data(_data.data().filter(function(d) {return d.time == _frame_number & d.dim == _dimension & d.selected===true}))
      .enter()
    .append("line")
      .attr("class", "dot_img")
      .attr("x1", function(d) { return d.b_x; })
      .attr("y1", function(d) { return d.b_y; })
      .attr("x2", function(d) { return d.d_x; })
      .attr("y2", function(d) { return d.d_y; })
      .style("stroke","cyan");
}

// runImageAnimation(_data, _frames, _dimension)
//   Input:
//     _data : list of all persistence points complete with selection status
//     _frames : a list of frame numbers (used for animation of images)
//     _dimension : an integer indicating dimension of persistence diagram of interest
//     _divs : the CSS selectors indicating the divs to paint animation
//     _height : the actual display height of the image
//     _width : the actual display width of the image
//   Effect:
//     Plays an animation of the images. Each image frame is annotated with
//     features selected with reverse selector brush.
function runImageAnimation (_data, _frames, _dimension, _divs, _height, _width) {

  var first_frame = _frames[0];
  var last_frame = _frames[_frames.length-1];
  var frame_index = 0;

  var playInterval = setInterval(function() {
    // Get current frame
    var current_frame = _frames[frame_index];

    // Draw current frame
    d3.select(_divs.divImg + ">img[frame_number='"+current_frame+"']")
      .attr("style","margin: 0; width: {0}px; height: {1}px;".format(_width, _height));

    // Clear previous frame
    var previous_frame = (frame_index == 0) ? last_frame : _frames[frame_index-1];
    if ( previous_frame != current_frame) {
      // Prevent erasing previous frame in one frame case
      d3.select(_divs.divImg + ">img[frame_number='"+previous_frame+"']")
        .attr("style","display: none;");
    }
    // Update slide number
    d3.select(_divs.divSlide + ">text")
      .text("Frame number: " + current_frame);

    // Update the drawings overlaid on the iamge
    annotateImageWithFeatures(_data, current_frame, d3.select(_divs.divImg+">svg"), _dimension);

    // Increment frame number
    ++ frame_index;

    // If frames exhausted, end animation
    if ( frame_index == _frames.length ) {
      clearInterval(playInterval); // setInterval returns and binds "playInterval" before this is encountered
    }
  }, 250);
};

// pointInPolygon (point, vs)
//   Input:
//     point : a point of the form [x, y]
//     vs : a list of vertices describing a polygon in the plane
//   Output:
//     Returns true if the point is within the polygon, false otherwise
//   Algorithm:
//     The algorithm computes the parity of the winding number of the curve
//     around the point by considering how many oriented segments of the polygon
//     intersect the horizontal rays from the point. 
function pointInPolygon (_point, _vs) {
  // from https://github.com/substack/point-in-polygon
  // ray-casting algorithm based on
  // http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
  var xi, xj, i, intersect, x = _point[0], y = _point[1], inside = false;
  // clever for loop to get each consecutive pair (i,j) circularly, i.e.
  //  (0,vs.length-1), (1, 0), (2, 1), ..., (vs.length-1, vs.length-2)
  for (var i = 0, j = _vs.length - 1; i < _vs.length; j = i++) { 
    xi = _vs[i][0], yi = _vs[i][1], xj = _vs[j][0], yj = _vs[j][1],
    intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
    if (intersect) inside = !inside;
  }
  return inside;
};

// plotPersistenceDiagram(_data, _frames, _dimension, _display_settings)
//   Input:
//     _data : a list of all points to be plotted in persistence diagram
//     _frames : a list of frame numbers (used for animation of images)
//     _dimension : an integer indicating dimension of persistence diagram of interest
//     _display_settings : information about display settings
//   Effect:
//     Creates a D3-interactive persistence diagram using the provided data.
//     The _display_settings input is an object containing display infomration, 
//     width, height, and functions getColor and getSize which are used to color 
//     the points.
function plotPersistenceDiagram( _data, _frames, _dimension, _display_settings ){

  // Pull out information from "_display_settings"
  var width = _display_settings.width;
  var height = _display_settings.height;
  var getColor = _display_settings.getColor;
  var getSize = _display_settings.getSize;
  var divs = _display_settings.divs;
  var pdsize = _display_settings.pdsize;
  var margin = 30; // Make room for the axis labels

  var x = d3.scale.linear()
      .range([0, pdsize]);

  var y = d3.scale.linear()
      .range([pdsize, 0]);

  var xAxis = d3.svg.axis()
      .scale(x)
      .orient("bottom");

  var yAxis = d3.svg.axis()
      .scale(y)
      .orient("left");

  var svg = d3.select(divs.divPD)
      .append("svg")
        .attr("width", pdsize + margin )
        .attr("height", pdsize + margin )
        .attr("class","svg")
        .append("g")
        .attr("transform", "translate(" + margin + ",0)")
        .style("font-size", "8pt");

  // Transforms placement of center of points onto plane (vs absolute coordinates)
  // x.domain(d3.extent(_data, function(d) { return d.birth; })).nice();
  // y.domain(d3.extent(_data, function(d) { return d.death; })).nice();
  x.domain([0,255]).nice(); // Fix coordinates to the entire range of grayscale values
  y.domain([0,255]).nice(); // Fix coordinates to the entire range of grayscale values

  // Style and build the persistence diagram
  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + pdsize + ")")
      .call(xAxis)
    .append("text")
      .attr("class", "label")
      .attr("x", pdsize)
      .attr("y", -6)
      .style("text-anchor", "end")
      .text("birth");
  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis)
    .append("text")
      .attr("class", "label")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text("death");

  // Plot the persistence plane
  svg.selectAll(".dot_")
     .data(_data.filter(function(d) {return d.dim == _dimension;}))
     .enter().append("circle")
       .attr("class", "dot_")
       .attr("cx", function(d) { return x(d.birth); })
       .attr("cy", function(d) { return y(d.death); })
       .style("fill", getColor)
       .attr("r", getSize);

  // // Title of plot
  // svg.append("text")
  //     .attr("class", "title")
  //     .attr("x", pdsize/2)
  //     .attr("y", 0 )
  //     .attr("text-anchor", "middle")
  //     .text("Dim " + _dimension);

  // Lasso functionality 
  var mysvg = d3.select(divs.divPD + " svg");
  var dot = d3.selectAll(divs.divPD + " .dot_");
  var g = mysvg.append("g");
  // D3 tool for creating SVG lines from coordinate data
  var line = d3.svg.line();

  // List of coordinates used by lasso drawing and dragging routines
  var coords = [];

  // Drawing lasso path
  var drawPath = function(terminator) {
    g.append("path").attr({
      d: line(coords)
    });
    if (terminator) {
      g.select(".terminator").remove();
      g.append("path").attr({
        class: "terminator",
        d: line([coords[0], coords[coords.length-1]])
      });
    }
  };

  // Code Executed when lasso drag operation begins
  var dragStart = function() {
    coords = [];
    g.selectAll("path").remove();
  };

  // Code exectuted as lasso dragging occurs
  var dragMove = function() {
    dot.classed("selected", false);
    coords.push(d3.mouse(this));
    dot.each(function(d, i) {
      d.selected = false;
      var p_x = parseFloat(d3.select(this).attr("cx"));
      var p_y = parseFloat(d3.select(this).attr("cy"));
      point = [p_x + margin, p_y ];
      if (pointInPolygon(point, coords)) {
        d3.select(this).classed("selected", true)
        d.selected = true;
      }
    });
    drawPath();
  };

  // Code executed when lasso dragging ends 
  var dragEnd = function() {
    drawPath(true);
    runImageAnimation(dot, _frames, _dimension, divs, height, width);
  };

  // D3 behavior encapsularing "dragStart", "dragMove", and "dragEnd":
  var drag = d3.behavior.drag()
            .on("dragstart", dragStart)
            .on("drag", dragMove)
            .on("dragend", dragEnd);

  // Associate D3 lasso dragging behavior with mysvg
  mysvg.call(drag);
}

// PersistenceExplorer(_imagefiles, _persistencefiles, _frames, _dimension, _divs, _imagesize, _maximagesize, _pdsize)
//   Inputs:
//     _imagefiles : a list of filenames containing image files
//     _persistencefiles : a list of filenames containing persistence diagram data
//     _frames : a list of frame numbers (frame indexing corresponds to _imagefiles and _persistencefiles indexing)
//     _dimension : an integer indicating dimension of persistence diagram of interest
//     _divs : the CSS selectors indicating the divs to draw app in
//     _imagesize : an array [width, height] holding the actual dimensions of the image
//     _maximagesize : an integer giving the maximum height/width of the displayed image in pixels 
//                     (can be larger than the acual image to magnify)
//     _pdsize : an integer giving the size of the persistence diagram, in pixels.
//   Effect:
//     Creates PersistenceExplorer application examining the frames of the
//     _imagefiles image files and _persistencefiles  data files which are
//     indicated in the list of frames _frames. Only persistence points 
//     of _dimension are considered. Image is scaled to size _maximagesize.
//     Persistence plane is also scaled to size _pdsize.
function PersistenceExplorer(_imagefiles, _persistencefiles, _frames, _dimension, _divs, _imagesize, _maximagesize, _pdsize) {

  // Display settings
  var divs = _divs || { divImg : "#divImg", divPD : "#divPD", divSlide : "#divSlide" };

  var imagewidth = _imagesize[0];
  var imageheight = _imagesize[1];
  var imagemaxsize = _maximagesize;
  var pdsize = _pdsize;

  var imagescale = imagemaxsize/Math.max(imagewidth, imageheight);

  var width = imagewidth*imagescale;
  var height = imageheight*imagescale;

  var brush = d3.svg.brush()
      .x(d3.scale.identity().domain([0, width]))  // brushwidth
      .y(d3.scale.identity().domain([0, height])); // brushheight
  var sizeUnselected = 2;
  var sizeSelected = 4;
  var colorInside = d3.scale.linear()
    .domain([_frames[0], _frames[_frames.length-1]])
    .range(["orange", "red"]);
  var colorOutside = d3.scale.linear()
    .domain([_frames[0], _frames[_frames.length-1]])
    .range(["gray", "blue"]);
  var getColor = function(_d) { return isFeatureSelected(_d,brush.extent()) ? colorInside(_d.time) : colorOutside(_d.time);};
  var getSize = function(_d) { return isFeatureSelected(_d,brush.extent()) ? sizeSelected : sizeUnselected; };

  var display_settings = {
    width : width,
    height : height,
    getColor : getColor,
    getSize : getSize,
    divs : divs,
    pdsize : pdsize,
    imagescale : imagescale
  };

  // Initialize the Image Data
  loadImages(_imagefiles, _frames, divs);

  // Display first frame
  d3.select(divs.divImg + " img[frame_number='"+_frames[0]+"']")
    .attr('style','width: {0}px; height: {1}px;'.format(width, height));

  // Create Frame Indicator
  d3.select(divs.divSlide)
    .append("text")
      .text("Frame number: " + _frames[0]);

  // Create SVG to draw image annotations in
  var svgImg = d3.select(divs.divImg)
      .append("svg")
        .attr("width", width )
        .attr("height", height )
        .attr("margin", 0)
        .attr("padding", 0);
        // .style("position", "absolute");

  // Setup callback for reverse selection
  brush.on("brushend", function () {
      // Unselect and unfill all persistence points in diagram
      alldots = d3.select(divs.divPD+">svg").selectAll("circle")
            .style("fill",null)
            .attr("r",sizeUnselected);
      // Redraw persistence points in diagram
      alldots.each(function(d,j) {
        d3.select(this).style("fill", getColor);
        d3.select(this).attr("r", getSize);
      });
    });

  // Create SVG group and associate brush to it
  svgImg.append("g")
    .attr("class", "brush")
    .call(brush);

  // Initialize the Persistence Diagram 
  loadPersistenceData(_persistencefiles, _frames, imagescale).then(function(data){
    plotPersistenceDiagram(data, _frames, _dimension, display_settings);
  });

};
