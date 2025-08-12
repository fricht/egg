#ifndef POLYGONER_H
#define POLYGONER_H

#include "utils.h"

typedef struct {
    double x;
    double y;
} Point;
typedef Point Vec; // alias

typedef struct {
    uint npts;
    Point* shape;
} Polygon;

Polygon* interactive_polygon_selector();

#endif
