console.log('runner loaded')
MathJax.Hub.Config({
    tex2jax: {inlineMath: [["$&","&$"]]}
  });
/*
	needed functions:
	something to handle a click to run action and set initial conditions
	something to evolve a given state forward a given period of time
	Xsomething to convert a reduced phase state to a configuration state
	Xsomething render a configuration state
		geometry notes:
			L=1 -> 100 px
			circles are 100px wide
			geometry is 400x400
			canvas is
*/
var goTime = function(){
	console.log('going')
	if(running){return}
	running = true
	var start = readForm()
	var state0 = new Vector(start)
	drawLoop(state0)
}

// drawLoop(new Vector([0,0,0,1.001,0.05, -.1,.2,0.0]))
var drawLoop = function(startVec){
    var prec = 10**-7
    var maxSteps = 1000
    var tF = .035
    var pathOut = utils.RK45BouncingPath(
        utils.reduced_dipole_equations, startVec, 0., tF,
        prec,maxSteps)
    var stateF = pathOut[pathOut.length-1][1]
    draw(stateF.values)
	setTimeout(drawLoop.bind(null, stateF), 35)
}

document.getElementById('runner').addEventListener('click', goTime)
running = false

var readForm = function(){
	var form = document.getElementsByClassName('controls')[0]
	console.log(form)
	var initVals = new Array(8)
	initVals[0] = 1*form['phi_d'].value
	initVals[1] = 1*form['phi_t'].value
	initVals[2] = 1*form['tht'].value
	initVals[3] = 1*form['r'].value
	initVals[0+4] = 1*form['pd'].value
	initVals[1+4] = 1*form['pt'].value
	initVals[2+4] = 1*form['ptht'].value
	initVals[3+4] = 1*form['pr'].value
	console.log(initVals)
	draw(initVals)
	return initVals
	// debugger
}
var dumby = [.15,0, 0.,1.5, 0,0,0,0]
var draw = function(state){
	console.log('draw attempted')
	var canvas = document.getElementById("render_area")
	var context = canvas.getContext('2d')
	var width = canvas.width
	var height = canvas.height
	// console.log(width, height)
	var configState = phaseConvert(state)
	var circrad = 50
	var rad = circrad*configState[3]
	var tht = configState[2]
	var x1 = rad*cos(tht+Math.PI) + width/2
	var y1 = height-(rad*sin(tht+Math.PI) + height/2)
	var x2 = rad*cos(tht) + width/2
	var y2 = height-(rad*sin(tht) + height/2)
	var phi1 = configState[0]
	var phi2 = configState[1]
	context.fillStyle = '#fff'
	context.fillRect(0,0, width, height)
	drawSphere(context, x1,y1, phi1)
	drawSphere(context, x2,y2, phi2)

	// console.log(canvas)
}

var drawSphere = function(context, x, y, phi){
	context.fillStyle = '#000'
	context.strokeStyle = "#000";
	var circrad = 50
	context.beginPath()
	// console.log('r0', x,y)
	context.arc(x,y, circrad, 0, Math.PI*2)
	context.stroke()
	context.beginPath()
	context.moveTo(x,y)
	var endx = x+cos(phi)*circrad*.75
	var endy = (y-sin(phi)*circrad*.75)
	// console.log('mu',endx,endy)
	context.lineTo(endx, endy)
	context.stroke()

}

var phaseConvert = function(phaseState){
	var config = new Array(8)
	config[0] = (phaseState[1]+phaseState[0])/2
	config[1] = (phaseState[1]-phaseState[0])/2
	config[2] = phaseState[2]
	config[3] = phaseState[3]
	config[0+4] = (phaseState[0+4]+phaseState[1+4])/10
	config[1+4] = (phaseState[0+4]-phaseState[1+4])/10
	config[2+4] = phaseState[6]/phaseState[3]**2
	config[3+4] = phaseState[7]
	return config
}

