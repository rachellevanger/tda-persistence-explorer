// PersistenceExplorer.js
// MIT LICENSE. Rachel Levanger 2016.
// Contributors: Rachel Levanger (author) 
//               Shaun Harker (contributor)

// Requires: d3

// loadImages(_imagefilenames, _frames)
//   Inputs:
//     _imagefilenames : a list of filenames pointing to image data
//     _frames : a list of frames (corresponding to _imagefilenames indexing) of interest
//   Effects:
//     Loads the images _imagefilenames[_frames[i]] for 0 <= i < _frames.length.
//     More specifically, it creates img tags with src files in _imagefilenames 
//     for the indicated _frames.
//     It expects there is a "div" named "imgContainer"
//     It creates an <img> tag inside this div with id "img_#" 
//     (e.g img_1, img_2,...) for each frame_number in _frames
//     and sets its "src" attribute to the filename _imagefilenames[frame_number] 
function loadImages(_imagefilenames, _frames) {
  _frames.forEach(function(frame_number){
      d3.select("#imgContainer")
          .append('img')
            .attr('src', _imagefilenames[frame_number])
            .attr('id', 'img_'+frame_number)
            .attr('class', 'img')
            .attr("style","display: none");
  });
}

// loadPersistenceData(_files, _frames)
//   Inputs:
//     _files : a list of filenames pointing to persistence data.
//              Each file is a .csv file which contains fields for 
//              "birth" "death" "b_x" "b_y" "d_x" "d_y" and "dim"
//     _frames : a list of frames (corresponding to _files indexing) of interest
//   Outputs:
//     Loads the files _files[_frames[i]] for 0 <= i < _frames.length, and returns
//     a promise-wrapped list of objects with "birth" "death" "time" and "selected" fields.
//     The "time" field gives the index of the file in which the data was loaded.
//     The "selected" field is initialized as "false"
function loadPersistenceData(_files, _frames) {
  function PromiseToLoadSingleFile(frame_number){
    return new Promise(function(resolve, reject){
        var img = new Image()
        d3.csv(_files[frame_number], function(error, dataset) {
          dataset.forEach(function(d) {
            d.birth = parseInt(d.birth);
            d.death = parseInt(d.death);
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
  // TODO: should b_x and friends be pre-parsed? Would this effect other code?
  var b_x = parseFloat(_d.b_x);
  var b_y = parseFloat(_d.b_y);
  var d_x = parseFloat(_d.d_x);
  var d_y = parseFloat(_d.d_y);
  return ( ( _extent[0][0] <= b_x ) && ( b_x <= (_extent[1][0]) ) && ( _extent[0][1] <= b_y ) && ( b_y <= _extent[1][1] )
    || ( _extent[0][0] <= d_x ) && ( d_x <= (_extent[1][0]) ) && ( _extent[0][1] <= d_y ) && ( d_y <= _extent[1][1] ) );
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
      .attr("id",function(d,i) {return "dot_" + i;}) // added
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
      .attr("id",function(d,i) {return "dot_" + i;}) // added
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
//   Effect:
//     Plays an animation of the images. Each image frame is annotated with
//     features selected with reverse selector brush.
function runImageAnimation (_data, _frames, _dimension) {

  var first_frame = _frames[0];
  var last_frame = _frames[_frames.length-1];
  var frame_index = 0;

  var playInterval = setInterval(function() {
    // Get current frame
    var current_frame = _frames[frame_index];

    // Draw current frame
    d3.select("#img_"+current_frame)
      .attr("style","");

    // Clear previous frame
    var previous_frame = (frame_index == 0) ? last_frame : _frames[frame_index-1];
    d3.select("#img_"+previous_frame)
      .attr("style","display: none");

    // Update slide number
    d3.select("#slideNo")
      .text("Sample point: " + current_frame);

    // Update the drawings overlaid on the iamge
    annotateImageWithFeatures(_data, current_frame, d3.select("#svgImg"), _dimension);

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
//     The _display_settings input is an object containing margin infomration, 
//     width, height, and functions getColor and getSize which are used to color 
//     the points.
function plotPersistenceDiagram( _data, _frames, _dimension, _display_settings ){

  // Pull out information from "_display_settings"
  var margin = _display_settings.margin;
  var width = _display_settings.width;
  var height = _display_settings.height;
  var getColor = _display_settings.getColor;
  var getSize = _display_settings.getSize;

  var x = d3.scale.linear()
      .range([0, width]);

  var y = d3.scale.linear()
      .range([height, 0]);

  var xAxis = d3.svg.axis()
      .scale(x)
      .orient("bottom");

  var yAxis = d3.svg.axis()
      .scale(y)
      .orient("left");

  var svg = d3.select("#divPD")
      .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .attr("class","svg")
        .attr("id","svgPD")
        .append("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // Transforms placement of center of points onto plane (vs absolute coordinates)
  x.domain(d3.extent(_data, function(d) { return d.birth; })).nice();
  y.domain(d3.extent(_data, function(d) { return d.death; })).nice();

  // Style and build the persistence diagram
  console.log("style and build persistence diagram");
  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
    .append("text")
      .attr("class", "label")
      .attr("x", width)
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

  console.log("plot PD");
  console.log(getColor);
  console.log(getSize);
  console.log(_display_settings);
  // Plot the persistence plane
  svg.selectAll(".dot_")
     .data(_data.filter(function(d) {return d.dim == _dimension;}))
     .enter().append("circle")
       .attr("id",function(d,i) {return "dot_" + i;})
       .attr("class", "dot_")
       .attr("cx", function(d) { return x(d.birth); })
       .attr("cy", function(d) { return y(d.death); })
       .style("fill", getColor)
       .attr("r", getSize);

  console.log("title PD");

  // Title of plot
  svg.append("text")
      .attr("class", "title")
      .attr("x", width/2)
      .attr("y", 0 - (margin.top / 2))
      .attr("text-anchor", "middle")
      .text("Dim " + _dimension);

  // Lasso functionality 
  console.log("build lasso");
  var mysvg = d3.select("#svgPD");
  var dot = d3.selectAll(".dot_");
  var g = mysvg.append("g")
               .attr("id","pointselector");
  console.log("hello?");
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
      g.select("#terminator").remove();
      g.append("path").attr({
        id: "terminator",
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
    //console.log(coords);
    dot.each(function(d, i) {
      d.selected = false;
      var p_x = parseFloat(d3.select(this).attr("cx"));
      var p_y = parseFloat(d3.select(this).attr("cy"));
      point = [p_x + margin.left, p_y + margin.top];
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
    runImageAnimation(dot, _frames, _dimension);
    //if (listGens) {
    //  listGenerators(dot);
    //}
  };

  // D3 behavior encapsularing "dragStart", "dragMove", and "dragEnd":
  var drag = d3.behavior.drag()
            .on("dragstart", dragStart)
            .on("drag", dragMove)
            .on("dragend", dragEnd);

  // Associate D3 lasso dragging behavior with mysvg
  mysvg.call(drag);
}

// PersistenceExplorer(_imagefiles, _persistencefiles, _frames, _dimension)
//   Inputs:
//     _imagefiles : a list of filenames containing image files
//     _persistencefiles : a list of filenames containing persistence diagram data
//     _frames : a list of frame numbers (frame indexing corresponds to _imagefiles and _persistencefiles indexing)
//     _dimension : an integer indicating dimension of persistence diagram of interest
//   Effect:
//     Creates PersistenceExplorer application examining the frames of the
//     _imagefiles image files and _persistencefiles  data files which are
//     indicated in the list of frames _frames. Only persistence points 
//     of _dimension are considered
function PersistenceExplorer(_imagefiles, _persistencefiles, _frames, _dimension) {

  // Display settings
  var margin = {top: 20, right: 20, bottom: 30, left: 40};
  var width = 421 - margin.left - margin.right;
  var height = 421 - margin.top - margin.bottom;
  var brush = d3.svg.brush()
      .x(d3.scale.identity().domain([0, 421]))  // brushwidth
      .y(d3.scale.identity().domain([0, 421])); // brushheight
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
    margin : margin,
    width : width,
    height : height,
    getColor : getColor,
    getSize : getSize
  };

  // Initialize the Image Data
  loadImages(_imagefiles, _frames);

  // Display first frame
  d3.select("#img_"+_frames[0])
    .attr('style','');

  // Create Frame Indicator
  var svgSlide = d3.select("#slide")
      .append("text")
        .attr("id","slideNo")
        .text("Sample point: " + _frames[0]);

  // Create SVG to draw image annotations in
  var svgImg = d3.select("#divImg")
      .append("svg")
        .attr("width", width + margin.left + margin.right + 20)
        .attr("height", height + margin.top + margin.bottom + 20)
        .attr("id","svgImg")

  // Setup callback for reverse selection
  brush.on("brushend", function () {
      // Unselect and unfill all persistence points in diagram
      alldots = d3.select("#svgPD").selectAll("circle")
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
  loadPersistenceData(_persistencefiles, _frames).then(function(data){
    plotPersistenceDiagram(data, _frames, _dimension, display_settings);
  });

};
