#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

typedef unsigned int uint;

typedef struct {
    double x;
    double y;
} Point;
typedef Point Vec; // alias

#define NPTS 6

Point PTS[NPTS] = {
    {0, 4},
    {0, 0},
    {3, 0},
    {3, 1},
    {1, 1},
    {1, 4}
};

double area(Point** points, uint length) {
    /*
    compute area delimited by points
    points : array of references of Point
    length : length of the array
    */
    double sum = 0;
    // my formula : signed area of triangles using successive determinant
    //                                                        (<=> signed area of parallelogram)
    /* for (uint i = 1; i < length - 1; i++) {
        sum += (points[i]->x - points[0]->x) * (points[i + 1]->y - points[i]->y)
            - (points[i + 1]->x - points[i]->x) * (points[i]->y - points[0]->y);
    } */
    // shoelace formula (smarter, more efficient, see wikipedia)
    for (uint i = 0; i < length; i++) {
        uint j = (i + 1) % length;
        sum += points[i]->x * points[j]->y - points[j]->x * points[i]->y;
    }
    return sum / 2;
}

bool is_point_active(Point* line_point, Vec* line_normal, Point* point) {
    // checks if point is on the normal side of the line
    Vec vector = {point->x - line_point->x, point->y - line_point->y};
    double scalar_prod = line_normal->x * vector.x + line_normal->y * vector.y;
    return scalar_prod >= 0;
}

Point line_segment_intersection(Point* pt, Vec* n, Point* a, Point* b) {
    // lets hope this will never need any refactor
    double ncX = n->x * pt->x;
    double ncY = n->y * pt->y;
    double naX = n->x * a->x;
    double naY = n->y * a->y;
    double nbX = n->x * b->x;
    double nbY = n->y * b->y;
    double t = (ncX + ncY - naX - naY) / (nbX - naX + nbY - naY);
    double x = a->x + t * (b->x - a->x);
    double y = a->y + t * (b->y - a->y);
    Point intersection = {x, y};
    return intersection;
}

int sorting_x(const void* p1, const void* p2) {
    const Point* a = *(const Point**)p1;
    const Point* b = *(const Point**)p2;
    return a->x < b->x;
}
int sorting_y(const void* p1, const void* p2) {
    const Point* a = *(const Point**)p1;
    const Point* b = *(const Point**)p2;
    return a->y < b->y;
}

double split_figure(
    Point* center, Vec* normal, // [pointer] [pointer] line definition
    Point* crossings, // [array] intersections holder (given, to stay in scope)
    Point** computed_shape, // [array [pointer]] active points (forming the shape)
    uint* shape_size // [pointer] nomber of points in computed_shape
) {
    uint crossings_counter = 0;
    *shape_size = 0;
    bool was_active = is_point_active(center, normal, &PTS[0]);
    if (was_active) {
        computed_shape[*shape_size] = &PTS[0];
        (*shape_size)++;
    }
    for (uint c = 1; c < NPTS; c++ /* haha what a joke XD */ ) {
        bool is_active = is_point_active(center, normal, &PTS[c]);
        if (was_active != is_active) {
            crossings[crossings_counter] = line_segment_intersection(center, normal, &PTS[c - 1], &PTS[c]);
            computed_shape[*shape_size] = &crossings[crossings_counter];
            (*shape_size)++;
            crossings_counter++;
        }
        if (is_active) {
            computed_shape[*shape_size] = &PTS[c];
            (*shape_size)++;
        }
        was_active = is_active;
    }
    bool is_active = is_point_active(center, normal, &PTS[0]);
    if (was_active != is_active) {
        crossings[crossings_counter] = line_segment_intersection(center, normal, &PTS[NPTS - 1], &PTS[0]);
        computed_shape[*shape_size] = &crossings[crossings_counter];
        (*shape_size)++;
        crossings_counter++;
    }
    assert(crossings_counter % 2 == 0);
    Point* all_crossings[crossings_counter];
    for (uint i = 0; i < crossings_counter; i++) {
        all_crossings[i] = &crossings[i];
    }
    if (fabs(normal->x) >= fabs(normal->y)) {
        qsort(all_crossings, crossings_counter, sizeof(Point*), sorting_x);
    } else {
        qsort(all_crossings, crossings_counter, sizeof(Point*), sorting_y);
    }
    double cross_section = 0;
    for (uint i = 0; i < crossings_counter; i += 2) {
        double dx = all_crossings[i]->x - all_crossings[i + 1]->x;
        double dy = all_crossings[i]->y - all_crossings[i + 1]->y;
        cross_section += sqrt(dx * dx + dy * dy);
    }
    return cross_section;
}

Point compute_c_given_n(
    Vec* n, double* learn_rate, // actual data, vector expected normalized
    Point* mean_point, double* half_area, // cached data
    Point** shape, Point* crossings, uint* shape_size // reuse to avoid re-alloc
) {
    Point c = {mean_point->x, mean_point->y};
    double last_cross_section = 0;
    double err = *half_area; // start with max error
    while (1 /* i feel like a real C dev by putting 1 instead of true */) {
        Point new_c = c; // TODO : compute new_c
        if (last_cross_section > 0) {
            double movment = - *learn_rate * err / last_cross_section;
            new_c.x += n->x * movment;
            new_c.y += n->y * movment;
        }
        last_cross_section = split_figure(&new_c, n, crossings, shape, shape_size);
        double new_area = area(shape, *shape_size);
        double new_err = *half_area - new_area;
        if (fabs(new_err) >= fabs(err) || err == 0) {
            break;
        }
        err = new_err;
        c = new_c;
    }
    return c;
}

int main() {
    uint worst_case_shape = (uint)(NPTS + NPTS / 2);
    uint worst_case_crossings = NPTS;
    Point* shape[worst_case_shape];
    Point crossings[worst_case_crossings];
    uint shape_size;
    Point center = {3, 0};
    Vec normal = {1, 0};
    double learn_rate = 0.001;
    // compute cache data
    double Sx = 0;
    double Sy = 0;
    Point* all_area[NPTS];
    for (uint i = 0; i < NPTS; i++) {
        Sx += PTS[i].x;
        Sy += PTS[i].y;
        all_area[i] = &PTS[i];
    }
    Point mean_point = {Sx / NPTS, Sy / NPTS};
    double half_area = area(all_area, NPTS) / 2;
    Point c = compute_c_given_n(&normal, &learn_rate, &mean_point, &half_area, shape, crossings, &shape_size);
    printf("c = {%lf, %lf}", c.x, c.y);
}
