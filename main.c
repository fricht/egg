#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#define PI 3.141592653689793

typedef unsigned int uint;

typedef struct {
    double x;
    double y;
} Point;
typedef Point Vec; // alias

typedef struct {
    uint npts;
    Point* shape;
} Polygon;

// L-shape
// #define NPTS 6
// Point PTS[NPTS] = {
//     {0, 7},
//     {0, 0},
//     {7, 0},
//     {7, 1},
//     {1, 1},
//     {1, 7}
// };

// E
#define NPTS 12
Point PTS[NPTS] = {
    {0, 9},
    {0, 0},
    {6, 0},
    {6, 1},
    {1, 1},
    {1, 4},
    {5, 4},
    {5, 5},
    {1, 5},
    {1, 8},
    {6, 8},
    {6, 9}
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

Point line_line_intersection(Point* c1, Vec* n1, Point* c2, Vec* n2) {
    // if it works, don't touch it (it doesn't)
    double q1 = n1->x * c1->x + n1->y * c1->y;
    double q2 = n2->x * c2->x + n2->y * c2->y;
    double denomX = n2->x * n1->y - n1->x * n2->y;
    double x = (n1->y * q2 - n2->y * q1) / denomX;
    double y = (n1->x * q2 - n2->x * q1) / -denomX;
    Point intersection = {x, y};
    return intersection;
}

int sorting_x(const void* p1, const void* p2) {
    const Point* a = *(const Point**)p1;
    const Point* b = *(const Point**)p2;
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    return 0;
}
int sorting_y(const void* p1, const void* p2) {
    const Point* a = *(const Point**)p1;
    const Point* b = *(const Point**)p2;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    return 0;
}

double split_figure(
    Point* center, Vec* normal, // [pointer] [pointer] line definition
    Polygon* polygon, // [pointer] the polygon (as in the name)
    Point* crossings, // [array] intersections holder (given, to stay in scope)
    Point** computed_shape, // [array [pointer]] active points (forming the shape)
    uint* shape_size // [pointer] nomber of points in computed_shape
) {
    uint crossings_counter = 0;
    *shape_size = 0;
    bool was_active = is_point_active(center, normal, &polygon->shape[0]);
    if (was_active) {
        computed_shape[*shape_size] = &polygon->shape[0];
        (*shape_size)++;
    }
    for (uint c = 1; c < polygon->npts; c++ /* haha what a joke XD */ ) {
        bool is_active = is_point_active(center, normal, &polygon->shape[c]);
        if (was_active != is_active) {
            crossings[crossings_counter] = line_segment_intersection(center, normal, &polygon->shape[c - 1], &polygon->shape[c]);
            computed_shape[*shape_size] = &crossings[crossings_counter];
            (*shape_size)++;
            crossings_counter++;
        }
        if (is_active) {
            computed_shape[*shape_size] = &polygon->shape[c];
            (*shape_size)++;
        }
        was_active = is_active;
    }
    bool is_active = is_point_active(center, normal, &polygon->shape[0]);
    if (was_active != is_active) {
        crossings[crossings_counter] = line_segment_intersection(center, normal, &polygon->shape[polygon->npts - 1], &polygon->shape[0]);
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
    Polygon* polygon, // polygon
    Vec* n, double* learn_rate, // actual data, vector expected normalized
    Point* mean_point, double* half_area, // cached data
    Point** shape, Point* crossings, uint* shape_size // reuse to avoid re-alloc
) {
    Point c = {mean_point->x, mean_point->y};
    double last_cross_section = 0;
    double err = *half_area; // start with max error
    while (1 /* i feel like a real C dev by putting 1 instead of true */) {
        Point new_c = c;
        if (last_cross_section > 0) {
            double movment = - *learn_rate * err / last_cross_section;
            new_c.x += n->x * movment;
            new_c.y += n->y * movment;
        }
        last_cross_section = split_figure(&new_c, n, polygon, crossings, shape, shape_size);
        double new_area = area(shape, *shape_size);
        double new_err = *half_area - new_area;
        if (fabs(new_err) >= fabs(err)) {
            break;
        }
        err = new_err;
        c = new_c;
    }
    return c;
}

void you_spin_me_right_round(Polygon* polygon, Point* Cs, Vec* Ns, uint samples, double learn_rate) {
    // compute cache data
    double Sx = 0;
    double Sy = 0;
    Point* all_area[polygon->npts];
    for (uint i = 0; i < polygon->npts; i++) {
        Sx += polygon->shape[i].x;
        Sy += polygon->shape[i].y;
        all_area[i] = &polygon->shape[i];
    }
    Point mean_point = {Sx / polygon->npts, Sy / polygon->npts};
    double half_area = area(all_area, polygon->npts) / 2;
    // buffers & things
    uint worst_case_shape = (uint)(polygon->npts + polygon->npts / 2);
    uint worst_case_crossings = polygon->npts;
    Point* shape[worst_case_shape];
    Point crossings[worst_case_crossings];
    uint shape_size;
    double delta_angle = PI / samples;
    double angle = 0;
    for (uint i = 0; i < samples; i++) {
        Vec normal = {cos(angle), sin(angle)};
        Point c = compute_c_given_n(polygon, &normal, &learn_rate, &mean_point, &half_area, shape, crossings, &shape_size);
        Cs[i] = c;
        Ns[i] = normal;
        angle += delta_angle;
    }
}

void compute_enveloppe(Point* Cs, Vec* Ns, uint samples, Point* enveloppe) {
    for (uint i = 0; i < samples; i++) {
        uint j = (i + 1) % samples;
        enveloppe[i] = line_line_intersection(&Cs[i], &Ns[i], &Cs[j], &Ns[j]);
    }
}

void create_svg(const char* filename, Polygon* polygon, Point* enveloppe, uint envlp_size, bool show_envlp_points, bool show_envlp_lines) {
    FILE* svg_file = fopen(filename, "w");
    if (!svg_file) {
        printf("fopen '%s' failed", filename);
        return;
    }
    // compute canvas boundaries
    double minX = polygon->shape[0].x;
    double maxX = polygon->shape[0].x;
    double minY = polygon->shape[0].y;
    double maxY = polygon->shape[0].y;
    for (uint i = 1; i < polygon->npts; i++) {
        minX = (polygon->shape[i].x < minX) ? polygon->shape[i].x : minX;
        maxX = (polygon->shape[i].x > maxX) ? polygon->shape[i].x : maxX;
        minY = (polygon->shape[i].y < minY) ? polygon->shape[i].y : minY;
        maxY = (polygon->shape[i].y > maxY) ? polygon->shape[i].y : maxY;
    }
    // 10% margin
    double marginX = 0.1 * (maxX - minX);
    double marginY = 0.1 * (maxY - minY);
    minX -= marginX;
    maxX += marginX;
    minY -= marginY;
    maxY += marginY;
    // headers
    fprintf(svg_file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(
        svg_file,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\" viewBox=\"%lf %lf %lf %lf\">\n",
        minX, minY, maxX - minX, maxY - minX
    );
    double stroke_width = 0.04; // TODO : dynamically compute that
    // draw shape
    for (uint i = 0; i < polygon->npts; i++) {
        uint j = (i + 1) % polygon->npts;
        fprintf(
            svg_file,
            "<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" stroke=\"black\" stroke-width=\"%lf\"/>\n",
            polygon->shape[i].x, polygon->shape[i].y, polygon->shape[j].x, polygon->shape[j].y, stroke_width
        );
    }
    // draw enveloppe
    for (uint i = 0; i < envlp_size; i++) {
        uint j = (i + 1) % envlp_size;
        if (show_envlp_points) {
            fprintf(
                svg_file,
                "<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" fill=\"purple\"/>",
                enveloppe[i].x, enveloppe[i].y, stroke_width
            );
        }
        if (show_envlp_lines) {
            fprintf(
                svg_file,
                "<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" stroke=\"red\" stroke-width=\"%lf\"/>\n",
                enveloppe[i].x, enveloppe[i].y, enveloppe[j].x, enveloppe[j].y, stroke_width / 2
            );
        }
    }
    // close everything
    fprintf(svg_file, "</svg>");
    fclose(svg_file);
}

int main() {
    uint samples = 100;
    Point Cs[samples];
    Vec Ns[samples];
    double learn_rate = 0.01;
    Polygon polygon = {NPTS, PTS};
    you_spin_me_right_round(&polygon, Cs, Ns, samples, learn_rate);
    Point enveloppe[samples];
    compute_enveloppe(Cs, Ns, samples, enveloppe);
    create_svg("out.svg", &polygon, enveloppe, samples, true, true);
}
