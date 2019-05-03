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
		var delta1 = kern1.mult(25/216).add(
			kern3.mult(1408/2565)).add(
			kern4.mult(2197/4104)).add(
			kern5.mult(-1/5)
			)
		var delta2 = kern1.mult(16/135).add(
			kern3.mult(6656/12825)).add(
			kern3.mult(28561/56430)).add(
			kern3.mult(-9/50)).add(
			kern3.mult(2/55)
			)
		return [delta1, delta2]
	}
	API.stepScaler = function(prec, step, deltaO4, deltaO5){
		var rescale = 10**6
		var numer = prec*step
		var diffVec = Vector.sub(deltaO5, deltaO4)
		var proposed
	  for(var i=0; i<diffVec.values.length; i++){
	    var diff = Math.abs(diffVec.values[i])
	    if(diff===0){
	    	proposed = 1
	    	continue
	    }
	    proposed = (numer/(2*diff))**.25
	    if(proposed < rescale){
	    	rescale = proposed
	    }
	  };
    if(proposed < rescale){
    	rescale = proposed
    }
    return rescale
	}
	return API
})()
