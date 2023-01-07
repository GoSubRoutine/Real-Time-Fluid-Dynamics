/**
 Real-Time Fluid Dynamics for Games by Jos Stam:
 https://DGP.Toronto.edu/public_user/stam/reality/Research/pub
 https://DGP.Toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf

 p5js: VrtXArt (@vrtx)
 https://GitHub.com/VrtXArt/cover3
 https://VrtXArt.GitHub.io/cover3
 
 mod: @GoToLoop (2023/Jan/06) [v2.1.5]
 https://GitHub.com/GoSubRoutine/Real-Time-Fluid-Dynamics
 https://GoSubRoutine.GitHub.io/Real-Time-Fluid-Dynamics

 https://Discourse.Processing.org/t/code-work-on-browser-too-slow-why/40455/3
*/

"use strict";

const
  FPS = 'FPS: ',
  ITERS = 20,

  N = 128,
  M = N + 1,
  NN = (N + 4) ** 2, // 17424
  NNN = NN << 1, // 34848

  SOURCE = 10,
  SOURCE2 = SOURCE >> 1, // 5

  DIFF = 1e-4, // .0001
  VISC = 1e-4,
  DT = .01,
  RSTEP = .2,
  NOISE_AMOUNT = .03,

  u = new Float32Array(NN),
  u_prev = u.slice(),

  v = new Float32Array(NN),
  v_prev = v.slice(),

  dens = new Float32Array(NN),
  dens_prev = dens.slice(),

  turb = new Float32Array(NNN),
  next_turb = turb.slice(),

  temp = new Float32Array(2);

var img, song, squareColor, black, white, fg, rpos;

function preload() {
  img = loadImage `capapng.png`;
}

function setup() {
  createCanvas(1280, 600);
  pixelDensity(1);

  describe `Fluids gradually fade away like smoke`;

  song = loadSound `Data.mp3`;

  squareColor = color `magenta`;
  fg = black = color(0);
  white = color(255);

  img.resize(width, height);

  initSim();
}

function draw() {
  document.title = FPS + nf(frameRate(), 2, 1);

  background(img);

  u_prev.set(u);
  v_prev.set(v);
  dens_prev.set(dens);

  if (mouseIsPressed) {
    add_density();
    add_velocity();
  }

  vel_step();
  add_noise();
  dens_step();

  drawDensity();

  squareColor.setAlpha(sin(millis() / 900) * 30 + 30);
  fill(squareColor).rect(0, 0, width, height);

  fill(fg).square(1000, 530, 50, 20);
}

function mousePressed() {
  if (mouseButton == CENTER)  return initSim();

  if (mouseX >= 1000 && mouseX <= 1050 && mouseY >= 530 && mouseY <= 560) {
    song.isPlaying()? song.stop() : song.play();
    fg = fg == white && black || white;
  }
}

function IX(i, j) {
  return (N + 2) * j + i;
}

function PX(x, y) {
  return width * y + x << 2;
}

function initSim() {
  rpos = 0;

  u.fill(0);
  v.fill(0);
  dens.fill(0);

  for (var i = 0; i < NNN; i += 2) {
    turb.set(polarBoxMullerTransform(temp), i);
    next_turb.set(polarBoxMullerTransform(temp), i);
  }
}

/**
 Basic implementation of the polar form of the Box-Muller transform

 Returns an array containing two Gaussian distributed random values
 with mean 0 and a standard deviation of 1
*/
function polarBoxMullerTransform(arr = new Float32Array(2)) {
  var x1, x2, w;

  do {
    x1 = random(-1, 1);
    x2 = random(-1, 1);
    w = x1*x1 + x2*x2;
  } while (w >= 1);

  w = sqrt(-2 * log(w) / w);

  arr[0] = x1 * w;
  arr[1] = x2 * w;

  return arr;
}

function add_density() {
  const nw = N / width * mouseX | 0, nh = N / height * mouseY | 0;

  dens[IX(nw, nh - 1)] += SOURCE2;
  dens[IX(nw - 1, nh)] += SOURCE2;
  dens[IX(nw + 1, nh)] += SOURCE2;
  dens[IX(nw + 9, nh)] += SOURCE;
  dens[IX(nw, nh + 1)] += SOURCE2;
}

function add_velocity() {
  const
    nw = N / width,
    nh = N / height,

    xv = nw * (mouseX - pmouseX),
    yv = nh * (mouseY - pmouseY),

    i = IX(nw * mouseX | 0, nh * mouseY | 0);

  u[i] += xv * 95 * (2 / (abs(xv) + 1));
  v[i] += yv * 35 * (2 / (abs(yv) + 1));
}

function set_bnd(opt, arr) {
  for (var i = 1; i <= N; ++i) {
    arr[IX(i, 0)] = opt == 2? -arr[IX(i, 1)] : arr[IX(i, 1)];
    arr[IX(0, i)] = opt == 1? -arr[IX(1, i)] : arr[IX(1, i)];
    arr[IX(M, i)] = opt == 1? -arr[IX(N, i)] : arr[IX(N, i)];
    arr[IX(i, M)] = opt == 2? -arr[IX(i, N)] : arr[IX(i, N)];
  }

  arr[IX(0, 0)] = .5 * (arr[IX(1, 0)] + arr[IX(0, 1)]);
  arr[IX(M, 0)] = .5 * (arr[IX(N, 0)] + arr[IX(M, 1)]);
  arr[IX(0, M)] = .5 * (arr[IX(1, M)] + arr[IX(0, N)]);
  arr[IX(M, M)] = .5 * (arr[IX(N, M)] + arr[IX(M, N)]);
}

function diffuse(opt, a, b, diff) {
  const d = DT * N*N * diff, d4 = 4*d + 1;

  for (var k = 0; k < ITERS; ++k) {
    for (var j = 1; j <= N; ++j) {
      const jj = IX(0, j), j0 = IX(0, j - 1), j1 = IX(0, j + 1);

      for (var i = 1; i <= N; ++i) {
        const ij = jj + i;

        a[ij] = (b[ij] + (
          a[j0 + i] +
          a[ij - 1] +
          a[ij + 1] +
          a[j1 + i]
        ) * d) / d4;
      }
    }

    set_bnd(opt, a);
  }
}

function advect(opt, a, b, u0, v0) {
  const dt0 = DT * N, lim = N + .5;

  for (var j = 1; j <= N; ++j) {
    const jj = IX(0, j);

    for (var i = 1; i <= N; ++i) {
      const
        ij = jj + i,

        x = constrain(i - dt0 * u0[ij], .5, lim),
        y = constrain(j - dt0 * v0[ij], .5, lim),

        i0 = ~~x, i1 = i0 + 1,
        j0 = ~~y, j1 = j0 + 1,

        s1 = x - i0, s0 = 1 - s1,
        t1 = y - j0, t0 = 1 - t1;

      a[ij] =
        s0 * (t0 * b[IX(i0, j0)] + t1 * b[IX(i0, j1)]) +
        s1 * (t0 * b[IX(i1, j0)] + t1 * b[IX(i1, j1)]);
    }
  }

  set_bnd(opt, a);
}

function project(u0, v0, p, div) {
  const h = 1 / N, h1 = -.5 * h, h2 = 2 * h;

  for (var j = 1; j <= N; ++j) {
    const jj = IX(0, j), j0 = IX(0, j - 1), j1 = IX(0, j + 1);

    for (var i = 1; i <= N; ++i) {
      const ij = jj + i;

      p[ij] = 0;

      div[ij] = h1 * (
        u0[ij + 1] - u0[ij - 1] +
        v0[j1 + i] - v0[j0 + i]
      );
    }
  }

  set_bnd(0, div);
  set_bnd(0, p);

  for (var k = 0; k < ITERS; ++k) {
    for (j = 1; j <= N; ++j) {
      const jj = IX(0, j), j0 = IX(0, j - 1), j1 = IX(0, j + 1);

      for (i = 1; i <= N; ++i) {
        const ij = jj + i;

        p[ij] = .25 * (
          div[ij] +
          p[j0 + i] +
          p[ij - 1] +
          p[ij + 1] +
          p[j1 + i]
        );
      }
    }

    set_bnd(0, p);
  }

  for (j = 1; j <= N; ++j) {
    const jj = IX(0, j), j0 = IX(0, j - 1), j1 = IX(0, j + 1);

    for (i = 1; i <= N; ++i) {
      const ij = jj + i;

      u0[ij] -= (p[ij + 1] - p[ij - 1]) / h2;
      v0[ij] -= (p[j1 + i] - p[j0 + i]) / h2;
    }
  }

  set_bnd(1, u0);
  set_bnd(2, v0);
}

function vel_step() {
  diffuse(1, u_prev, u, VISC);
  diffuse(2, v_prev, v, VISC);
  project(u_prev, v_prev, u, v);
  advect(1, u, u_prev, u_prev, v_prev);
  advect(2, v, v_prev, u_prev, v_prev);
  project(u, v, u_prev, v_prev);
}

function add_noise() {
  if ((rpos += RSTEP) >= 1)  rpos = 0;

  const rpos1 = 1 - rpos;

  for (var y = 1; y <= N; ++y) {
    const yy = IX(0, y), y0 = IX(0, y - 1), y1 = IX(0, y + 1);

    for (var x = 1; x <= N; ++x) {
      const ii = yy + x, i = ii << 1, j = i + 1;

      rpos || refill(i, j);

      const
        hg = abs(dens[ii - 1] - dens[ii + 1]),
        vg = abs(dens[y0 + x] - dens[y1 + x]),

        un = hg * (turb[i] * rpos1 + next_turb[i] * rpos),
        vn = vg * (turb[j] * rpos1 + next_turb[j] * rpos);

      u[ii] += un * NOISE_AMOUNT * (2 / (abs(un) + 1));
      v[ii] += vn * NOISE_AMOUNT * (2 / (abs(vn) + 1));
    }
  }
}

function refill(i, j = i + 1) {
  turb[i] = next_turb[i];
  turb[j] = next_turb[j];

  next_turb.set(polarBoxMullerTransform(temp), i);
}

function dens_step() {
  diffuse(0, dens_prev, dens, DIFF);
  advect(0, dens, dens_prev, u, v);
}

function drawDensity() {
  const NW = N / width, NH = N / height;

  loadPixels();

  for (var y = 0; y < height; ++y) {
    const
      dy = NH * y,
      dyc = ceil(dy),
      dyf = ~~dy,
      ddy = dy - dyf,
      ddy1 = 1 - ddy;

    for (var x = 0; x < width; ++x) {
      const
        dx = NW * x,
        dxc = ceil(dx),
        dxf = ~~dx,
        ddx = dx - dxf,
        ddx1 = 1 - ddx,

        df =
          ddy1 * (ddx1 * dens[IX(dxf, dyf)] + ddx * dens[IX(dxc, dyf)]) +
          ddy  * (ddx1 * dens[IX(dxf, dyc)] + ddx * dens[IX(dxc, dyc)]),

        di = constrain(df * 255 | 0, 0, 255),

        xy = PX(x, y);

      pixels[xy] = pixels[xy] * (1 - df) + df*di;
      pixels[xy + 1] = pixels[xy + 2] = di;
      pixels[xy + 3] = 255;
    }
  }

  updatePixels();
}
