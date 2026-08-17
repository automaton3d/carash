/*
 * carash.c
 *
 * CaRaSh-style local spherical shell propagation test.
 *
 * Goal: expand a wavefront from a central source in a 3-D cubic lattice so
 * that, at each tick t, exactly the cells with integer Euclidean radius r = t
 * become active.  No sqrt, no isqrt, no general multiplication in the inner
 * rule; only additions and shifts are used.
 *
 * Each cell carries:
 *   ax, ay, az  signed offset from the source (propagated locally)
 *   r2          squared Euclidean distance (r2 = ax^2 + ay^2 + az^2)
 *   r           integer radius (floor(sqrt(r2)))
 *   next_sq     (r+1)^2, updated by adding 2*r+1 (shift + 1)
 *   tick        the tick at which this cell becomes active ( == r )
 *
 * The source emits at tick 0.  A cell that becomes active at tick t tries
 * every "outward" 6-neighbour direction (the direction that increases the
 * absolute value of the corresponding coordinate, or both directions when the
 * coordinate is zero).  The neighbour's r2 is updated by
 *     r2_new = r2 + 2*|coord| + 1
 * which is a shift and an addition, not a multiplication.  The neighbour's
 * new radius is either r or r+1, decided by comparing r2_new with next_sq.
 *
 * Because adjacent cells on the same integer-radius shell can activate each
 * other within the same tick (tangential propagation), each main tick is
 * iterated to a fixed point before moving to t+1.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifndef L
#define L 17
#endif

#define CENTER (L / 2)
#define INF    0x7FFFFFFF
#define MAXT   (4 * L)

typedef struct {
    int ax, ay, az;          /* signed offset from source */
    unsigned int r2;       /* squared Euclidean distance */
    int r;                 /* integer radius */
    unsigned int next_sq;  /* (r+1)^2, kept incrementally */
    int tick;              /* activation tick */
    int fired;             /* last main tick in which this cell fired */
} Cell;

static Cell grid[L][L][L];
static int  exact[L][L][L];

/* Reference radius: floor(sqrt(dx^2+dy^2+dz^2)).
 * This is allowed to use sqrt because it is only the ground truth, not the
 * algorithm under test. */
static int ref_radius(int dx, int dy, int dz)
{
    double d = sqrt((double)dx * dx + (double)dy * dy + (double)dz * dz);
    return (int)d;
}

static void init_exact(void)
{
    int x, y, z;
    for (x = 0; x < L; ++x)
        for (y = 0; y < L; ++y)
            for (z = 0; z < L; ++z)
                exact[x][y][z] = ref_radius(x - CENTER, y - CENTER, z - CENTER);
}

static void init_grid(void)
{
    int x, y, z;
    for (x = 0; x < L; ++x)
        for (y = 0; y < L; ++y)
            for (z = 0; z < L; ++z) {
                Cell *c = &grid[x][y][z];
                c->ax = c->ay = c->az = 0;
                c->r2 = 0;
                c->r = -1;
                c->next_sq = 0;
                c->tick = INF;
                c->fired = -1;
            }

    Cell *src = &grid[CENTER][CENTER][CENTER];
    src->ax = src->ay = src->az = 0;
    src->r2 = 0;
    src->r = 0;
    src->next_sq = 1;   /* (0+1)^2 = 1 */
    src->tick = 0;
    src->fired = -1;
}

/* Try to update one outward neighbour.  sx, sy, sz are one of the 6 cardinal
 * directions.  Only directions that increase |ax|, |ay| or |az| are used. */
static void try_propagate(int x, int y, int z, int sx, int sy, int sz)
{
    int nx = x + sx;
    int ny = y + sy;
    int nz = z + sz;
    if (nx < 0 || nx >= L || ny < 0 || ny >= L || nz < 0 || nz >= L)
        return;

    Cell *c = &grid[x][y][z];
    Cell *n = &grid[nx][ny][nz];

    /* Choose the coordinate that changes. */
    int coord = 0;
    int new_ax = c->ax, new_ay = c->ay, new_az = c->az;
    if (sx != 0) {
        coord = c->ax;
        new_ax = c->ax + sx;
    } else if (sy != 0) {
        coord = c->ay;
        new_ay = c->ay + sy;
    } else {
        coord = c->az;
        new_az = c->az + sz;
    }

    int abs_coord = coord < 0 ? -coord : coord;
    /* diff = 2*|coord| + 1  (a shift plus one) */
    unsigned int diff = ((unsigned int)abs_coord << 1) + 1u;
    unsigned int new_r2 = c->r2 + diff;

    int new_r;
    unsigned int new_next_sq;
    if (new_r2 < c->next_sq) {
        new_r = c->r;
        new_next_sq = c->next_sq;
    } else {
        new_r = c->r + 1;
        /* next_sq = (r+2)^2 = (r+1)^2 + 2*(r+1) + 1 */
        new_next_sq = c->next_sq + ((unsigned int)new_r << 1) + 1u;
    }

    int new_tick = c->tick + (new_r - c->r);  /* 0 or 1 */

    int better = 0;
    if (new_tick < n->tick)
        better = 1;
    else if (new_tick == n->tick && new_r2 < n->r2)
        better = 1;

    if (!better)
        return;

    n->ax = new_ax;
    n->ay = new_ay;
    n->az = new_az;
    n->r2 = new_r2;
    n->r = new_r;
    n->next_sq = new_next_sq;
    n->tick = new_tick;
}

/* Fire one cell: propagate to every outward 6-neighbour. */
static void fire(int x, int y, int z, int t)
{
    Cell *c = &grid[x][y][z];
    if (c->fired == t)
        return;
    c->fired = t;

    /* x direction */
    if (c->ax >= 0 && x + 1 < L) try_propagate(x, y, z, +1, 0, 0);
    if (c->ax <= 0 && x - 1 >= 0) try_propagate(x, y, z, -1, 0, 0);
    /* y direction */
    if (c->ay >= 0 && y + 1 < L) try_propagate(x, y, z, 0, +1, 0);
    if (c->ay <= 0 && y - 1 >= 0) try_propagate(x, y, z, 0, -1, 0);
    /* z direction */
    if (c->az >= 0 && z + 1 < L) try_propagate(x, y, z, 0, 0, +1);
    if (c->az <= 0 && z - 1 >= 0) try_propagate(x, y, z, 0, 0, -1);
}

/* One main tick: propagate to fixed point (handles tangential same-tick
 * activation inside the shell). */
static int tick(int t)
{
    int changed;
    do {
        changed = 0;
        int x, y, z;
        for (x = 0; x < L; ++x)
            for (y = 0; y < L; ++y)
                for (z = 0; z < L; ++z) {
                    Cell *c = &grid[x][y][z];
                    if (c->tick == t && c->fired != t) {
                        fire(x, y, z, t);
                        changed = 1;
                    }
                }
    } while (changed);
    return 0;
}

static int verify(void)
{
    int errors = 0;
    int x, y, z;
    for (x = 0; x < L; ++x)
        for (y = 0; y < L; ++y)
            for (z = 0; z < L; ++z) {
                Cell *c = &grid[x][y][z];
                if (c->tick == INF)
                    continue;
                if (c->r != exact[x][y][z]) {
                    printf("  MISMATCH at (%d,%d,%d): algorithm r=%d, exact r=%d, tick=%d, r2=%u\n",
                           x, y, z, c->r, exact[x][y][z], c->tick, c->r2);
                    errors++;
                }
                if (c->tick != exact[x][y][z]) {
                    printf("  MISMATCH at (%d,%d,%d): algorithm tick=%d, exact tick=%d\n",
                           x, y, z, c->tick, exact[x][y][z]);
                    errors++;
                }
            }
    return errors;
}

static void print_slice(int z)
{
    printf("z = %d slice (ticks == exact radius):\n", z);
    int x, y;
    for (y = 0; y < L; ++y) {
        for (x = 0; x < L; ++x) {
            Cell *c = &grid[x][y][z];
            if (c->tick == INF)
                printf("  .");
            else
                printf(" %2d", c->tick);
        }
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    int max_t = 0;
    (void)argc; (void)argv;

    init_exact();
    init_grid();

    printf("CaRaSh local spherical shell test (L=%d, center=%d)\n", L, CENTER);
    printf("Running ticks until the grid is filled...\n");

    int t;
    for (t = 0; t < MAXT; ++t) {
        tick(t);

        /* Stop early if the outer sphere has been reached everywhere. */
        int all_filled = 1;
        int x, y, z;
        for (x = 0; x < L; ++x)
            for (y = 0; y < L; ++y)
                for (z = 0; z < L; ++z)
                    if (grid[x][y][z].tick == INF && exact[x][y][z] <= t + 1)
                        all_filled = 0;

        if (all_filled) {
            max_t = t;
            break;
        }
    }

    printf("Stopped at tick %d\n", max_t);

    int errors = verify();
    printf("Verification: %d mismatches\n", errors);

    /* Per-tick shell population. */
    printf("\nShell population per tick (active cells at each tick):\n");
    int pop[MAXT];
    memset(pop, 0, sizeof(pop));
    int x, y, z;
    for (x = 0; x < L; ++x)
        for (y = 0; y < L; ++y)
            for (z = 0; z < L; ++z)
                if (grid[x][y][z].tick >= 0 && grid[x][y][z].tick < MAXT)
                    pop[grid[x][y][z].tick]++;

    for (t = 0; t <= max_t; ++t)
        printf("  tick %2d: %4d cells\n", t, pop[t]);

    print_slice(CENTER);

    return errors ? 1 : 0;
}
