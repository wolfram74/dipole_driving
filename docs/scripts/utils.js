cos = Math.cos
sin = Math.sin
utils = (function(){
	var API = {}
	API.rk45Steps = function(gradFunc, t, state, stepSize){
		var kern1 = gradFunc(t, state).mult(stepSize);
		var kern2 = gradFunc(
			t+stepSize/4,
			state.add(kern1.mult(1/4))
			).mult(stepSize)
		var kern3 = gradFunc(
			t+3*stepSize / 8,
			state.add(kern1.mult(3/32)).add(kern2.mult(9/32))
			).mult(stepSize)
		var kern4 = gradFunc(
			t+12*stepSize/13,
			state.add(
				kern1.mult(1932/2197)
				).add(
				kern2.mult(-7200/2197)
				).add(
				kern3.mult(7296/2197)
				)
			).mult(stepSize)
		var kern5 = gradFunc(
			t+stepSize/1,
			state.add(
				kern1.mult(439/216)
				).add(
				kern2.mult(-8/1)
				).add(
				kern3.mult(3680/513)
				).add(
				kern4.mult(-845/4104)
				)
			).mult(stepSize)
		var kern6 = gradFunc(
			t+stepSize/2,
			state.add(
				kern1.mult(-8/27)
				).add(
				kern2.mult(2/1)
				).add(
				kern3.mult(-3544/2565)
				).add(
				kern4.mult(1859/4104)
				).add(
				kern5.mult(-11/40)
				)
			).mult(stepSize)
		// console.log(kern1.values, kern2.values, kern3.values, kern4.values, kern5.values, kern6.values)
		var delta1 = kern1.mult(25/216).add(
			kern3.mult(1408/2565)).add(
			kern4.mult(2197/4104)).add(
			kern5.mult(-1/5)
			)
		var delta2 = kern1.mult(16/135).add(
			kern3.mult(6656/12825)).add(
			kern4.mult(28561/56430)).add(
			kern5.mult(-9/50)).add(
			kern6.mult(2/55)
			)
		return [delta1, delta2]
	}

	API.RK45SimplePath = function(
		gradFunc, state, startTime, endTime,
		precision, maxSteps
		){
		var stepSize = 10**-2
		var path = [[startTime, state]]
		var running = true
		var lastLoop = false
		var timeLeft = endTime - startTime
		while(running){
			lastTime = path[path.length-1][0]
			lastState = path[path.length-1][1]
			if(lastLoop){
				running=false
			}
			if(path.length > maxSteps){
				running = false
				console.log('did not finish, time left', timeLeft)
			}

	    var steps = API.rk45Steps(gradFunc, lastTime, lastState, stepSize)
	    var stepO4 = steps[0]
	    var stepO5 = steps[1]
	    var stepAdjust = API.stepScaler(precision, stepSize, stepO5, stepO4)
	    // console.log(stepAdjust, lastTime)
	    if(stepAdjust < .9){
	    	stepSize *= stepAdjust
	    	continue
	    }
	    newState = lastState.add(stepO5)
	    newTime = lastTime+stepSize
	    path.push([newTime, newState ])

	    if(stepAdjust > 1.25){
	    	stepAdjust = 1.25
	    }
	    stepSize *= stepAdjust

	    timeLeft = endTime - newTime
	    if( timeLeft < stepSize){
	    	console.log('finishing early, steps taken=', path.length)
	    	stepSize = timeLeft
	    	lastLoop = true
	    }
		}
		return path
	}

	API.RK45BouncingPath = function(
		gradFunc, state, startTime, endTime,
		precision, maxSteps
		){
		var stepSize = 10**-2
		var path = [[startTime, state]]
		var running = true
		var lastLoop = false
		var timeLeft = endTime - startTime
		var boundary = 1.
		var margin = 10**-4
		while(running){
			lastTime = path[path.length-1][0]
			lastState = path[path.length-1][1]
			if(lastLoop){
				running=false
			}
			if(path.length > maxSteps){
				running = false
				console.log('did not finish, time left', timeLeft)
			}

	    var steps = API.rk45Steps(gradFunc, lastTime, lastState, stepSize)
	    var stepO4 = steps[0]
	    var stepO5 = steps[1]
	    // console.log(steps)
	    var stepAdjust = API.stepScaler(precision, stepSize, stepO5, stepO4)
	    // console.log(stepAdjust, lastTime)
	    if(stepAdjust < .9){
	    	stepSize *= stepAdjust
	    	continue
	    }
	    newState = lastState.add(stepO5)
	    newTime = lastTime+stepSize

	    if(newState.values[3] < (boundary+margin) && newState.values[7] <0){
	    	newState.values[7] *= -1
	    }

	    path.push([newTime, newState ])

	    if(stepAdjust > 1.25){
	    	stepAdjust = 1.25
	    }
	    stepSize *= stepAdjust

	    if(newState.values[7] < 0){
	    	tGuess = (boundary+margin-newState.values[3])/(2*newState.values[7])
	    	if( tGuess < stepSize){
	    		stepSize = tGuess
	    	}
	    }

	    timeLeft = endTime - newTime
	    if( timeLeft < stepSize){
	    	// console.log('finishing early, steps taken=', path.length)
	    	stepSize = timeLeft
	    	lastLoop = true
	    }
		}
		return path
	}

	API.stepScaler = function(prec, step, deltaO4, deltaO5){
		var rescale = 10**6
		var numer = prec*step
		var diffVec = Vector.sub(deltaO5, deltaO4)
		// console.log('diff', diffVec.values)
		var proposed
	  for(var i=0; i<diffVec.values.length; i++){
	    var diff = Math.abs(diffVec.values[i])
	    // console.log(diff, diff==0)
	    if(diff===0){
	    	proposed = 1
	    	continue
	    }
	    proposed = (numer/(2*diff))**.25
	    // console.log(proposed)
	    if(proposed < rescale){
	    	rescale = proposed
	    }
	  };
    return rescale
	}

	API.reduced_dipole_equations = function(t, state){
		var deltas = new Array(state.values.length)
    //velocities
    deltas[0] = 20.* state.values[4]
    deltas[1] = 20.* state.values[5]
    deltas[2] = 2. * state.values[6]*(state.values[3]**(-2.))
    if(state.values[3]==1. && state.values[7]<0){
    	deltas[3]=0
    }else{deltas[3] = 2. * state.values[7]}
    //forces
    deltas[4] = -sin(state.values[0])/(12.*state.values[3]**3)
    deltas[5] = -sin(state.values[1] - 2.*state.values[2])/(4.*state.values[3]**3)
    deltas[6] = 2.*sin(state.values[1] - 2.*state.values[2])/(4.*state.values[3]**3)
    if( state.values[3]==1. && state.values[7]<0){
    	deltas[7] = 0
    }else{
	    deltas[7] = (
	        2.*state.values[6]**2/state.values[3]**3
	        -(
	            cos(state.values[0]) + 3.*cos(state.values[1] - 2.*state.values[2])
	        )/(4.*state.values[3]**4)
	    )
    }
    return new Vector(deltas)
	}

	API.totalEnergy = function(state){
		// console.log(API.PE(state), API.KE(state))
		return API.PE(state) + API.KE(state)
	}

	API.totalL = function(state){
    return 2.*state.values[5]+state.values[6]
	}

	API.PE = function(state){
    return -(
        cos(state.values[0]) + 3*cos(state.values[1] - 2.*state.values[2])
        )/(12.*state.values[3]**3)
	}

	API.KE = function(state){
		// console.log((2*state.values[7]**2)/2)
    return (
        2*state.values[7]**2
        +2*state.values[6]**2/state.values[3]**2
        +20*(state.values[4]**2+state.values[5]**2)
        )/2
	}


	return API
})()
