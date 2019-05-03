var Vector = function(array){
  this.values = array
}

Vector.prototype.add = function(otherVec){
  if(this.values.length != otherVec.values.length){
    throw new Error("Dimension mismatch on vector addition")
  };
  for(var i=0; i<this.values.length; i++){
    this.values[i] += otherVec.values[i]
  };
  return this
};

Vector.prototype.mult= function(scalar){
  for(var i=0; i<this.values.length; i++){
    this.values[i] *= scalar
  };
  return this
}

Vector.prototype.sub= function(otherVec){
  // untested
  return this.add( otherVec.mult(-1) )
}

Vector.prototype.div= function(scalar){
  // untested
  this.mult(1/scalar)
  return this
}

Vector.prototype.dot = function(otherVec){
  if(this.values.length != otherVec.values.length){
    throw new Error("Dimension mismatch on vector addition")
  };
  var sum = 0
  for(var i=0; i<this.values.length; i++){
    sum += this.values[i] * otherVec.values[i]
  };
  return sum
}

Vector.prototype.mag = function(){
  // untested
  return Math.pow(this.dot(this), 0.5)
}

Vector.prototype.setMag = function(scalar){
  // untested
  var magnitude = this.mag()
  return this.mult(scalar/magnitude)
};

Vector.prototype.normalize = function(){
  // untested
  var magnitude = this.mag()
  return this.div(magnitude)
}

Vector.prototype.limit = function(scalar){
  var magnitude = this.mag();
  if(magnitude > scalar){ return this.setMag(scalar)};
  return this
};

Vector.prototype.modulus = function(scalar){
  // fackin' joova scripts
  for(var i=0; i<this.values.length; i++){
    this.values[i] = ((this.values[i] % scalar)+scalar)%scalar
  };
  return this
}

//class functions

Vector.mult = function(initVector, scalar){
  var values = []
  for(var i=0; i<initVector.values.length; i++){
    values.push(initVector.values[i] * scalar)
  };
  return new Vector(values)
};

Vector.copy = function(initVector){
  return this.mult(initVector, 1)
}

Vector.add = function(vect1, vect2){
  var values = []
  for(var i=0; i<vect1.values.length; i++){
    values.push(vect1.values[i] + vect2.values[i])
  };
  return new Vector(values)
}

Vector.div = function(initVector, scalar){
  return this.mult(initVector, 1/scalar)
};

Vector.sub = function(vect1, vect2){
  return this.add(vect1, this.mult(vect2, -1))
};
