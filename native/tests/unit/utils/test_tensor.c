#include "utils/tensor.h"

#include <math.h>
#include <stddef.h>
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

static void test_zeroes_only_requested_contiguous_components(void)
{
    double storage[8] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    size_t index;

    rp_tensor_zero(&storage[2], 4);

    CHECK(storage[0] == 1.0);
    CHECK(storage[1] == 2.0);
    for (index = 2; index < 6; ++index) {
        CHECK(storage[index] == 0.0);
    }
    CHECK(storage[6] == 7.0);
    CHECK(storage[7] == 8.0);
}

static void test_zero_count_does_not_access_storage(void)
{
    double value = 7.0;

    rp_tensor_zero(&value, 0);
    rp_tensor_zero(NULL, 0);
    CHECK(value == 7.0);
}

static void test_finiteness_for_arbitrary_contiguous_regions(void)
{
    const double finite_values[4] = {-1.0, 0.0, 2.5, 1e300};
    const double nonfinite_values[4] = {1.0, NAN, 2.0, INFINITY};

    CHECK(rp_tensor_is_finite(finite_values, 4));
    CHECK(rp_tensor_is_finite(nonfinite_values, 1));
    CHECK(!rp_tensor_is_finite(nonfinite_values, 2));
    CHECK(!rp_tensor_is_finite(&nonfinite_values[2], 2));
    CHECK(rp_tensor_is_finite(NULL, 0));
}

int main(void)
{
    test_zeroes_only_requested_contiguous_components();
    test_zero_count_does_not_access_storage();
    test_finiteness_for_arbitrary_contiguous_regions();

    if (failures != 0) {
        (void)fprintf(stderr, "%d tensor utility assertion(s) failed\n", failures);
        return 1;
    }
    (void)puts("Tensor utility tests passed");
    return 0;
}
