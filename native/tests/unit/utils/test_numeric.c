#include "utils/numeric.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

static int failures = 0;

#define CHECK(condition)                                                        \
    do {                                                                        \
        if (!(condition)) {                                                     \
            (void)fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__,     \
                #condition);                                                    \
            ++failures;                                                         \
        }                                                                       \
    } while (0)

int main(void)
{
    CHECK(rp_effectively_zero(0.0, 0.0));
    CHECK(rp_effectively_zero(64.0 * DBL_EPSILON, 1.0));
    CHECK(rp_effectively_zero(-64.0 * DBL_EPSILON, 1.0));
    CHECK(!rp_effectively_zero(65.0 * DBL_EPSILON, 1.0));
    CHECK(!rp_effectively_zero(1.0, 0.0));
    CHECK(rp_bl_polar_axis_singular(64.0 * DBL_EPSILON));
    CHECK(!rp_bl_polar_axis_singular(65.0 * DBL_EPSILON));
    CHECK(rp_bl_polar_axis_singular(acos(-1.0) - 64.0 * DBL_EPSILON));
    CHECK(!rp_bl_polar_axis_singular(acos(-1.0) - 128.0 * DBL_EPSILON));
    CHECK(rp_bl_polar_axis_singular(-0.1));
    CHECK(rp_bl_polar_axis_singular(acos(-1.0) + 0.1));

    if (failures != 0) {
        (void)fprintf(stderr, "%d numeric utility assertion(s) failed\n", failures);
        return 1;
    }
    (void)puts("Numeric utility tests passed");
    return 0;
}
