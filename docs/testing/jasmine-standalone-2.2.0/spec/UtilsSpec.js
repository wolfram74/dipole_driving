describe("Meta Utilities", function() {
  it("should run tests", function(){
    expect(true).toEqual(true)
  })
});

function base_SHO(t, state){
  // state is a Vector object
  deltas = new Array(state.values.length)
  deltas[0] = state.values[1]
  deltas[1] = -state.values[0]
  return new Vector(deltas)
}


describe('utils behavior', function(){
  it('rk45 steps make two deltas', function(){
    var stateI = new Vector([0, 1])
    var step = 10**-4
    var prec = 10**-6
    var steps = utils.rk45Steps(base_SHO, 0, stateI, step)
    var stepO4 = steps[0]
    var stepO5 = steps[1]
    console.log(stepO4.values)
    console.log(stepO5.values)
    expect(stepO4.values.length).toEqual(stepO5.values.length)
    var adjust1 = utils.stepScaler(prec, step, stepO4, stepO5)
    step *= adjust1
    steps = utils.rk45Steps(base_SHO, 0, stateI, step)
    stepO4 = steps[0]
    stepO5 = steps[1]
    var adjust2 = utils.stepScaler(prec, step, stepO4, stepO5)
    console.log(adjust1, adjust2)
    // python version got (2788.318719217494, 0.5039921097456969)
    expect(Math.abs(Math.log(adjust2))).toBeLessThan(Math.abs(Math.log(adjust1)))
  })
  xit('rk45SimplePath works for a SHO', function(){
    var stateI = new Vector([0, 1])
    var prec = 10**-6
    var maxSteps = 10
    var tF = 2*Math.PI
    var pathOut = utils.RK45SimplePath(
      base_SHO, stateI, 0, tF,
      prec, maxSteps
      )
    var stateF = pathOut[pathOut.length-1]
    expect(Math.abs(stateF[0]-tF)).toBeLessThan(10**-5)
    expect(Math.abs(stateF[1].values[0])).toBeLessThan(10**-5)
  })
})
