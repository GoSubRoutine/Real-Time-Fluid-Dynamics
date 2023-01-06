// Source: Real-Time Fluid Dynamics for Games by Jos Stam - http://www.intpowertechcorp.com/GDC03.pdf
let song;
var overButton = false;
var myboolean = (backwards = false);
let value1 = 0;

var N = 128;
var size;
var u = [];
var v = [];
var u_prev = [];
var v_prev = [];
var dens = [];
var dens_prev = [];
var source = 10;
var diff = 0.0001;
var visc = 0.0001;
var dt = 0.01;
var turb = [];
var next_turb = [];
var rpos = 0;
var rstep = 0.2;
var noise_amount = 0.03;

let img;

function preload() {
  img = loadImage('capapng.png');
}

function IX( i = 0, j = 0 ) {
  return i + (N + 2) * j;
}

function PX( x = 0, y = 0 ) {
  return (x + width * y) * 4;
}

function setup() {
  createCanvas(1280, 600);
  pixelDensity(1);
  //background(0);
  //stroke(255,0,250);
  initSim();
  
  song = loadSound("Data.mp3");
}

function initSim() {
  size = (N + 4) * (N + 4);
  for (var i = 0; i < size; i++) {
    u[i] = 0.0;
    v[i] = 0.0;
    dens[i] = 0.0;
    turb[i] = polarBoxMullerTransform();
    next_turb[i] = polarBoxMullerTransform();
  }
}

function draw() {
  clear();
  background(0);
  
  document.title = 'FPS: ' + nf(frameRate(), 2, 1);
  
  dens_prev = dens.slice();
  u_prev = u.slice();
  v_prev = v.slice();
  add_density();
  add_velocity();
  vel_step();
  add_noise();
  dens_step();
  image(img, 0, 0);
  img.resize(1280, 600);
  //drawVelocity();
  drawDensity();

  squareColor = color(255, 0, 255);
  squareColor.setAlpha(30 + 30 * sin(millis() / 900));
  fill(squareColor);
  rect(0, 0, width, height );
  describe('a square with gradually changing opacity on a gray background');
  
  
    //buttons draws
  fill(value1);
  rect(1000, 530, 50, 50, 20);
}

  

function add_density() {
  if (mouseIsPressed) {

    dens[IX(int( (N / width) * mouseX+9 ), int( (N / height) * mouseY ))] += source;
    dens[IX(int( (N / width) * mouseX - 1), int( (N / height) * mouseY ))] += source / 2;
    dens[IX(int( (N / width) * mouseX + 1), int( (N / height) * mouseY ))] += source / 2;
    dens[IX(int( (N / width) * mouseX), int( (N / height) * mouseY - 1))] += source / 2;
    dens[IX(int( (N / width) * mouseX), int( (N / height) * mouseY + 1))] += source / 2;
  }
}

function add_velocity() {
  var i;
  if (mouseIsPressed) {
    i = IX(int( (N / width) * mouseX ), int( (N / height) * mouseY ));
    var xv = (N / width) * (mouseX - pmouseX);
    var yv = (N / height) * (mouseY - pmouseY);
    u[i] += xv * (2 / (abs(xv) + 1)) * 95;
    v[i] += yv * (2 / (abs(yv) + 1)) * 35;
  }
}

function add_noise() {

  var refill = false;
  rpos += rstep;
  if (rpos >= 1) {
    refill = true;
    rpos = 0;
  }
  for (var x=1; x<=N; x++) {
    for (var y=1; y<=N; y++) {
      var i = IX(x, y);
      if (refill) {
        turb[i] = next_turb[i];
        next_turb[i] = polarBoxMullerTransform();
      }
      var hg = abs( dens[IX(x-1, y)] - dens[IX(x+1, y)] );
      var vg = abs( dens[IX(x, y-1)] - dens[IX(x, y+1)] );
      var un = (turb[i][0] * (1.0 - rpos) + next_turb[i][0] * rpos) * hg;
      var vn = (turb[i][1] * (1.0 - rpos) + next_turb[i][1] * rpos) * vg;
      u[i] += un * (2 / (abs(un) + 1)) * noise_amount;
      v[i] += vn * (2 / (abs(vn) + 1)) * noise_amount;
    }
  }
}

function diffuse( b, x, x0, diff0 ) {
  var a = dt * diff0 * N * N;
  for ( var k=0; k<20; k++ ) {
    for ( var i=1; i<=N; i++ ) {
      for ( var j=1; j<=N; j++ ) {
        x[IX(i, j)] = ( x0[IX(i, j)] + a * (x[IX(i-1, j)] + x[IX(i+1, j)] + x[IX(i, j-1)] + x[IX(i, j+1)]) ) / (1 + 4 * a);
      }
    }
    set_bnd ( b, x );
  }
}

function advect( b, d, d0, u0, v0 ) {
  var i0, j0, i1, j1;
  var x, y, s0, t0, s1, t1, dt0;
  dt0 = dt * N;
  for ( var i=1; i<=N; i++ ) {
    for ( var j=1; j<=N; j++ ) {
      x = i - dt0 * u0[IX(i, j)];
      y = j - dt0 * v0[IX(i, j)];
      if (x < 0.5) x = 0.5;
      if (x > N + 0.5) x = N + 0.5;
      i0 = int(x);
      i1 = i0+1;
      if (y < 0.5) y = 0.5;
      if (y > N + 0.5) y = N + 0.5;
      j0 = int(y);
      j1 = j0 + 1;
      s1 = x - i0;
      s0 = 1 - s1;
      t1 = y - j0;
      t0 = 1 - t1;
      d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
    }
  }
  set_bnd ( b, d );
}

function dens_step() {
  //SWAP ( x0, x );
  diffuse( 0, dens_prev, dens, diff );
  //SWAP ( x0, x );
  advect( 0, dens, dens_prev, u, v );
}

function vel_step() {
  //SWAP( u0, u );
  diffuse( 1, u_prev, u, visc );
  //SWAP( v0, v );
  diffuse( 2, v_prev, v, visc );
  project( u_prev, v_prev, u, v );
  //SWAP( u0, u );
  //SWAP( v0, v );
  advect( 1, u, u_prev, u_prev, v_prev );
  advect ( 2, v, v_prev, u_prev, v_prev );
  project( u, v, u_prev, v_prev );
}

function project ( u0, v0, p, div ) {
  var h = 1.0 / N;
  for ( var i=1; i<=N; i++ ) {
    for ( var j=1; j<=N; j++ ) {
      div[IX(i, j)] = -0.5 * h * ( u0[IX(i+1, j)] - u0[IX(i-1, j)] + v0[IX(i, j+1)] - v0[IX(i, j-1)] );
      p[IX(i, j)] = 0;
    }
  }
  set_bnd ( 0, div );
  set_bnd ( 0, p );
  for ( var k=0; k<20; k++ ) {
    for ( i=1; i<=N; i++ ) {
      for ( j=1; j<=N; j++ ) {
        p[IX(i, j)] = ( div[IX(i, j)] + p[IX(i-1, j)] + p[IX(i+1, j)] + p[IX(i, j-1)] + p[IX(i, j+1)] ) * (1/4);
      }
    }
    set_bnd ( 0, p );
  }
  for ( i=1; i<=N; i++ ) {
    for ( j=1; j<=N; j++ ) {
      u0[IX(i, j)] -= 0.5 * ( p[IX(i+1, j)] - p[IX(i-1, j)] ) / h;
      v0[IX(i, j)] -= 0.5 * ( p[IX(i, j+1)] - p[IX(i, j-1)] ) / h;
    }
  }
  set_bnd ( 1, u0 );
  set_bnd ( 2, v0 );
}

function set_bnd ( b, x ) {
  for ( var i=1; i<=N; i++ ) {
    x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
    x[IX(N+1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
    x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
    x[IX(i, N+1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
  }
  x[IX(0, 0 )] = 0.5 * ( x[IX(1, 0 )] + x[IX(0, 1)] );
  x[IX(0, N+1)] = 0.5 * ( x[IX(1, N+1)] + x[IX(0, N )] );
  x[IX(N+1, 0 )] = 0.5 * ( x[IX(N, 0 )] + x[IX(N+1, 1)] );
  x[IX(N+1, N+1)] = 0.5 * ( x[IX(N, N+1)] + x[IX(N+1, N )] );
}


function drawDensity() {

  var dx, dy, ddx, ddy;
  var df, di;
  loadPixels();

  for (var x = 0; x < width; x++) {
    for (var y = 0; y < height; y++) {
      dx = (N / width) * x;
      ddx = dx - int(dx);
      dy = (N / height) * y;
      ddy = dy - int(dy);
      df = (dens[IX(floor(dx), floor(dy))] * (1.0 - ddx) + dens[IX(ceil(dx), floor(dy))] * ddx) * (1.0 - ddy) + (dens[IX(floor(dx), ceil(dy))] * (1.0 - ddx) + dens[IX(ceil(dx), ceil(dy))] * ddx) * ddy;
      di = int(df * 255);
      if (di < 0) di = 0;
      if (di > 255) di = 255;
      pixels[PX(x, y)] = pixels[PX(x, y)] * (1-df) + di*df;
      pixels[PX(x, y) + 1] = di;
      pixels[PX(x, y) + 2] = di;
      pixels[PX(x, y) + 3] = 255;
    }
  }
  updatePixels();
  fill(0, 140, 160, 100);
  noStroke();
  //rect(0,height * .6, width, height);
}


function keyPressed() {
  //save("img_" + month() + '-' + day() + '_' + hour() + '-' + minute() + '-' + second() + ".jpg");
}

function drawVelocity() {
  var sx = width / N;
  var sy = height / N;
  for (var x = 1; x <= N; x++) {
    for (var y = 1; y <= N; y++) {
      var i = IX(x, y);
      /*var b = 0;
       var r = int( u[i] * v[i] * 255);
       if (r < 0) {
       b = -r;
       r = 0;
       }
       if (r > 255) r = 255;
       if (b > 255) b = 255;
       stroke(r, 0, b);*/
      line(int((x - 0.5) * sx), int((y - 0.5) * sy), int((x - 0.5) * sx + u[i] * 50), int((y - 0.5) * sy + v[i] * 50));
    }
  }
}

//Basic implementation of the polar form of the Box-Muller transform
//Returns an array containing two gaussian distributed random values with mean 0 and a standard deviation of 1
function polarBoxMullerTransform() {
  var x1, x2, w, y1, y2;
  do {
    x1 = 2.0 * random() - 1.0;
    x2 = 2.0 * random() - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );

  w = Math.sqrt( (-2.0 * Math.log( w ) ) / w );
  y1 = x1 * w;
  y2 = x2 * w;
  return [y1, y2];
}

function mouseReleased() {
  loop();

  if (mouseX >= 1000 && mouseX <= 1050 && (mouseY >= 530) & (mouseY <= 560)) {
    //song.play();

    if (song.isPlaying()) {
      // .isPlaying() returns a boolean
      song.stop();
    } else {
      song.play();
    }

    if (value1 == 255) {
      value1 = 0;
    } else {
      value1 = 255;
    }
    //print("btn1 hit.");
    
  }
}
