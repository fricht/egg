#include "polygoner.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// L-shape
#define L_NPTS 6
Point L_PTS[L_NPTS] = {
    {0, 7},
    {0, 0},
    {7, 0},
    {7, 1},
    {1, 1},
    {1, 7}
};
Polygon L_polygon = {L_NPTS, L_PTS};

// E
#define E_NPTS 12
Point E_PTS[E_NPTS] = {
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
Polygon E_polygon = {L_NPTS, L_PTS};

Polygon* const_polygon_selector() {
    char choice;
    do {
        printf(
            "\nWhich polygon do you want ?\n"
            " 0 - L\n"
            " 1 - E\n"
            " |-> "
        );
        choice = getchar();
        while (getchar() != '\n');
    } while (choice != '0' && choice != '1');
    switch (choice) {
        case '0':
            return &L_polygon;
        case '1':
            return &E_polygon;
    }
    printf("Error, switch skipped. defaulting to 'L_polygon'\n");
    return &L_polygon;
}

Polygon* star_polygon_creator() {
    int branches;
    printf(
        "\nHow many branches do you want ? (>= 2)\n"
        " |-> "
    );
    if (scanf("%d", &branches) != 1) {
        printf("NAN Error. int expected\n");
        exit(1);
    }
    if (branches < 2) {
        printf("sorry, at least 2 branches are required\n");
        exit(1);
    }
    double branch_size;
    printf(
        "\nWhat size ?\n"
        " |-> "
    );
    if (scanf("%lf", &branch_size) != 1) {
        printf("NAN Error. double expected\n");
        exit(1);
    }
    // this part of the code will malloc but not free
    // its not a problem currently since when the data
    // isn't used anymore the program finishes, but
    // be careful re-using that code elsewhere
    uint npts = 3 * branches;
    Point* star_shape = malloc(npts * sizeof(Point)); // Warning !!!
    Point base_pts[branches];
    double delta_angle = 2 * PI / branches;
    double angle = 0;
    for (uint i = 0; i < branches; i++) {
        Point pt = {cos(angle), sin(angle)};
        base_pts[i] = pt;
        angle += delta_angle;
    }
    for (uint i = 0; i < branches; i++) {
        uint j = (i + 1) % branches;
        star_shape[3 * i] = base_pts[i];
        double ext_angle = (2 * i + 1) * PI / branches;
        Vec extender = {branch_size * cos(ext_angle), branch_size * sin(ext_angle)};
        Point ext1 = {base_pts[i].x + extender.x, base_pts[i].y + extender.y};
        Point ext2 = {base_pts[j].x + extender.x, base_pts[j].y + extender.y};
        star_shape[3 * i + 1] = ext1;
        star_shape[3 * i + 2] = ext2;
    }
    Polygon* star_polygon = malloc(sizeof(Polygon)); // Warning !!!
    star_polygon->npts = npts;
    star_polygon->shape = star_shape;
    return star_polygon;
}

Polygon* interactive_polygon_selector() {
    char choice;
    do {
        printf(
            "\nWhat kind of polygon do you want ?\n"
            " 0 - const\n"
            " 1 - star\n"
            " 2 - random\n"
            " |-> "
        );
        choice = getchar();
        while (getchar() != '\n');
    } while (choice != '0' && choice != '1' && choice != '2');
    switch (choice) {
        case '0':
            return const_polygon_selector();
        case '1':
            return star_polygon_creator();
    }
    printf("Invalid number. defaulting to 'L_polygon'\n");
    return &L_polygon;
}
