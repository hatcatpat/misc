#include <cairo/cairo.h>

#include <cglm/cglm.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265358979323846
#define TAU (2.0 * PI)
#define SQRT_2 1.4142135623730951

#define MIN(a, b) (a) < (b) ? (a) : (b)
#define MAX(a, b) (a) > (b) ? (a) : (b)

typedef unsigned int uint;

float f_rand(float lo, float hi) {
  return ((float)rand() / (float)RAND_MAX) * (hi - lo) + lo;
}
float norm_rand() { return f_rand(0.0, 1.0); }
float bi_rand() { return f_rand(-1.0, 1.0); }

float bi(float norm) { return norm * 2.0 - 1.0; }
float norm(float bi) { return (bi + 1.0) * 0.5; }

float norm_scale(float a, float lo, float hi) { return a * (hi - lo) + lo; }
float bi_scale(float a, float lo, float hi) {
  return (a + 1.0) * 0.5 * (hi - lo) + lo;
}

// noise taken from p5js :)
#define PERLIN_YWRAPB 4
#define PERLIN_YWRAP (1 << PERLIN_YWRAPB)
#define PERLIN_ZWRAPB 8
#define PERLIN_ZWRAP (1 << PERLIN_ZWRAPB)
#define PERLIN_SIZE 4095
float perlin_octaves = 16;      // default to medium smooth
float perlin_amp_falloff = 0.5; // 50% reduction/octave
float perlin[PERLIN_SIZE];
float scaled_cos(float theta) { return 0.5 * (1.0 - cos(theta * PI)); }
float noise(float x, float y, float z) {
  if (x < 0)
    x = -x;
  if (y < 0)
    y = -y;
  if (z < 0)
    z = -z;

  int xi = floor(x), yi = floor(y), zi = floor(z);
  float xf = x - xi;
  float yf = y - yi;
  float zf = z - zi;
  float rxf, ryf;

  float r = 0;
  float ampl = 0.5;

  float n1, n2, n3;

  for (int o = 0; o < perlin_octaves; o++) {
    int of = xi + (yi << PERLIN_YWRAPB) + (zi << PERLIN_ZWRAPB);

    rxf = scaled_cos(xf);
    ryf = scaled_cos(yf);

    n1 = perlin[of & PERLIN_SIZE];
    n1 += rxf * (perlin[(of + 1) & PERLIN_SIZE] - n1);
    n2 = perlin[(of + PERLIN_YWRAP) & PERLIN_SIZE];
    n2 += rxf * (perlin[(of + PERLIN_YWRAP + 1) & PERLIN_SIZE] - n2);
    n1 += ryf * (n2 - n1);

    of += PERLIN_ZWRAP;
    n2 = perlin[of & PERLIN_SIZE];
    n2 += rxf * (perlin[(of + 1) & PERLIN_SIZE] - n2);
    n3 = perlin[(of + PERLIN_YWRAP) & PERLIN_SIZE];
    n3 += rxf * (perlin[(of + PERLIN_YWRAP + 1) & PERLIN_SIZE] - n3);
    n2 += ryf * (n3 - n2);

    n1 += scaled_cos(zf) * (n2 - n1);

    r += n1 * ampl;
    ampl *= perlin_amp_falloff;
    xi <<= 1;
    xf *= 2;
    yi <<= 1;
    yf *= 2;
    zi <<= 1;
    zf *= 2;

    if (xf >= 1.0) {
      xi++;
      xf--;
    }
    if (yf >= 1.0) {
      yi++;
      yf--;
    }
    if (zf >= 1.0) {
      zi++;
      zf--;
    }
  }

  return r;
}

void cairo_set_source_hsv(cairo_t *cr, float h, float s, float v) {
  float c = v * s;
  float x = c * (1.0 - fabs(fmod(h / 60.0, 2.0) - 1.0));
  float m = v - c;

  int n = (int)floor(h / 60.0);
  float r, g, b;
  switch (n) {
  case 0:
    r = c, g = x, b = 0.0;
    break;

  case 1:
    r = x, g = c, b = 0.0;
    break;

  case 2:
    r = 0.0, g = c, b = x;
    break;

  case 3:
    r = 0.0, g = x, b = c;
    break;

  case 4:
    r = x, g = 0.0, b = c;
    break;

  case 5:
    r = c, g = 0.0, b = x;
    break;
  }

  r += m, g += m, b += m;

  cairo_set_source_rgb(cr, r, g, b);
}

//
// config //
//

// - size
#define WIDTH 800
#define HEIGHT 800
#define MARGIN 0.9

// - drawing

//#define DRAW_GRID
#define RADIUS 0.3

//#define USE_ROUNDED
//#define LINE_WIDTH 0.005
#define LINE_WIDTH 0.0011
//#define LINE_WIDTH 0.01

#define DOT_SIZE 0.0009

//#define USE_DASH
#define DASH_SIZE 1.0

// - flow
#define GRID_SIZE 64
#define MAX_STEPS 1024
#define FLOW_NOISE_SCALE 4.0

//#define USE_GRID_POS

#define USE_CHANCE
#define CHANCE 0.4

//
// variables //
//

// - cairo
cairo_t *cr;
cairo_surface_t *surface;

// - flow
typedef struct {
  vec2 pos;
  char taken;
} cell;

cell grid[GRID_SIZE][GRID_SIZE];
float inv_grid_size, radius, diam;

typedef struct {
  vec2 pos;
  char active;
} walker;

walker walkers[GRID_SIZE][GRID_SIZE];

// - misc
typedef uint color[3];
color colors[] = {
    {0xff, 0xff, 0xff}, {0xff, 0, 0},
    //{0, 0xff, 0x00},
};
size_t num_colors = sizeof(colors) / sizeof(color);

//
//
//
void global_setup() {
  srand(time(NULL));

  surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, WIDTH, HEIGHT);
  cr = cairo_create(surface);

#ifdef USE_DASH
  double dashes[2] = {DASH_SIZE * LINE_WIDTH, DASH_SIZE * LINE_WIDTH};
  cairo_set_dash(cr, dashes, 2, 0.0);
#endif

#ifdef USE_ROUNDED
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
#endif

  cairo_scale(cr, WIDTH, HEIGHT);

  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_rectangle(cr, 0.0, 0.0, 1.0, 1.0);
  cairo_fill(cr);

  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_set_line_width(cr, LINE_WIDTH);

  cairo_scale(cr, MARGIN, MARGIN);
  float inv_margin = (1.0 - MARGIN) * 0.5;
  cairo_translate(cr, inv_margin, inv_margin);
}

void cleanup() {
  cairo_surface_write_to_png(surface, "image.png");

  cairo_destroy(cr);
  cairo_surface_destroy(surface);
}

void setup() {
  inv_grid_size = 1.0 / GRID_SIZE;
  radius = RADIUS * 0.5 * inv_grid_size;
  diam = radius * 2.0;

  cairo_translate(cr, inv_grid_size * 0.5, inv_grid_size * 0.5);

  for (int i = 0; i < GRID_SIZE; ++i) {
    for (int j = 0; j < GRID_SIZE; ++j) {
      {
        cell *c = &grid[i][j];
        c->pos[0] = (float)i * inv_grid_size;
        c->pos[1] = (float)j * inv_grid_size;
        c->taken = 0;
      }

      {
        walker *w = &walkers[i][j];

#ifdef USE_GRID_POS
        w->pos[0] = grid[i][j].pos[0];
        w->pos[1] = grid[i][j].pos[1];
#else
        w->pos[0] = f_rand(0.0, 1.0);
        w->pos[1] = f_rand(0.0, 1.0);
#endif

#ifdef USE_CHANCE
        w->active = f_rand(0.0, 1.0) < CHANCE;
#else
        w->active = 1;
#endif
      }
    }
  }
}

uint nearest_grid(float a) { return (uint)floor(a * GRID_SIZE); }

uint is_coord_valid(float a) {
  if (a < 0.0)
    return 0;
  else if (a > 1.0)
    return 0;

  return 1;
}

float flow(vec2 p) {
  vec2 p2;
  glm_vec2_copy(p, p2);
  glm_vec2_scale(p2, FLOW_NOISE_SCALE, p2);

  float d = f_rand(0.0, 0.01);

  float n = noise(p2[0], p2[1], 0);
  float f = tan((d + n) * TAU * 0.4);
  f = norm(f);

  return f;
}

void set_color(uint index) {
  color *col = &colors[index];

  float out[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    out[i] = (float)(*col)[i] / 255;

  cairo_set_source_rgb(cr, out[0], out[1], out[2]);
}

void loop() {
  vec2 p, q;
  walker *w;

  for (int i = 0; i < GRID_SIZE; ++i) {
    for (int j = 0; j < GRID_SIZE; ++j) {
      w = &walkers[i][j];

      if (!w->active)
        continue;

      float nsc = 8;
      float ip = i * inv_grid_size * nsc;
      float jp = j * inv_grid_size * nsc;
      float c = noise(0.0, ip, jp);
      c *= num_colors;
      set_color((uint)floor(c));

      float th = noise(100.0, ip, jp);
      th = norm_scale(th, 1.0, 4.0);
      cairo_set_line_width(cr, LINE_WIDTH * th);

      glm_vec2_copy(w->pos, p);
      cairo_move_to(cr, p[0], p[1]);

      for (int s = 1; s < MAX_STEPS; ++s) {
        printf("(%i, %i, %i)\n", i, j, s);

        glm_vec2_copy(w->pos, p);
        cairo_move_to(cr, p[0], p[1]);

        float f = flow(p);

        float c = noise(10.0, p[0] * 4 * FLOW_NOISE_SCALE,
                        p[1] * 4 * FLOW_NOISE_SCALE);
        if (c < 0.4) {
          p[0] = MIN(p[0] + f_rand(0.0, 8.0) * inv_grid_size, 1.0);
        } else if (c > 0.6) {
          p[1] = MIN(p[1] + f_rand(0.0, 8.0) * inv_grid_size, 1.0);
        }

        q[0] = p[0] + sin(f * TAU) * inv_grid_size;
        if (!is_coord_valid(q[0]))
          break;

        q[1] = p[1] + cos(f * TAU) * inv_grid_size;
        if (!is_coord_valid(q[1]))
          break;

        uint gx = nearest_grid(q[0]);
        uint gy = nearest_grid(q[1]);
        if (grid[gx][gy].taken)
          break;

        grid[gx][gy].taken = 1;
        glm_vec2_copy(q, w->pos);
        if (noise(0, i, j) < 0.5) {
          cairo_arc(cr, w->pos[0], w->pos[1], DOT_SIZE * th, 0, TAU);
          cairo_fill(cr);
        } else {
          cairo_line_to(cr, w->pos[0], w->pos[1]);
          cairo_stroke(cr);
        }
      }
    }
  }
}

void draw_grid() {
  for (int i = 0; i < GRID_SIZE; ++i) {
    for (int j = 0; j < GRID_SIZE; ++j) {
      cell *c = &grid[i][j];
      float f = flow(c->pos);

      float bf = norm_scale(f / TAU, 0.0, 1.0);
      cairo_set_source_rgb(cr, bf, bf, bf);

      cairo_arc(cr, c->pos[0], c->pos[1], radius, 0.0, TAU);
      cairo_fill(cr);

      cairo_move_to(cr, c->pos[0], c->pos[1]);
      cairo_line_to(cr, c->pos[0] + sin(f) * inv_grid_size * 0.8,
                    c->pos[1] + cos(f) * inv_grid_size * 0.8);
      cairo_stroke(cr);
    }
  }
}

int main(void) {
  global_setup();

  for (int i = 0; i < PERLIN_SIZE; ++i)
    perlin[i] = f_rand(0.0, 1.0);

  setup();
#ifdef DRAW_GRID
  draw_grid();
#endif
  loop();

  cleanup();
  return 0;
}
