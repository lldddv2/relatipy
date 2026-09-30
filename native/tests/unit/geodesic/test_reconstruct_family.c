#include "geodesic/solution/reconstruct.h"

#include <assert.h>
#include <math.h>
#include <string.h>

int main(void)
{
    double states[3][8] = {
        {0.0, 8.0, 1.1, 20.0, 1.2, -0.01, 0.002, 0.02},
        {0.0, 8.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.02},
        {0.0, 8.0, 1.1, 20.0, 1.2, -0.01, 0.002, 0.02}
    };
    double original[3][8];
    double full[3][29];
    double cartesian[3][7];
    double selected[3][7];
    double elements[3][6];
    rp_kerr_status statuses[3];
    size_t i;
    size_t family;
    memcpy(original, states, sizeof(states));
    assert(rp_solution_reconstruct_canonical_batch(0.5, &states[0][0], 3U,
        &full[0][0], &cartesian[0][0], statuses)
        == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    for (family = RP_SOLUTION_FAMILY_CARTESIAN;
        family <= RP_SOLUTION_FAMILY_SPHERICAL; ++family) {
        assert(rp_solution_reconstruct_canonical_family_batch(0.5,
            &states[0][0], 3U, (rp_solution_family)family, &selected[0][0],
            statuses) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
        assert(statuses[0] == RP_KERR_STATUS_OK);
        assert(statuses[1] == RP_KERR_STATUS_COORDINATE_SINGULARITY);
        assert(statuses[2] == RP_KERR_STATUS_OK);
        for (i = 0U; i < 7U; ++i) {
            assert(selected[1][i] == 0.0);
            assert(selected[0][i] == selected[2][i]);
            if (family == RP_SOLUTION_FAMILY_CARTESIAN) {
                assert(selected[0][i] == cartesian[0][i]);
            } else if (i == 0U) {
                assert(selected[0][i] == full[0][RP_SOL_T]);
            } else if (i < 4U) {
                assert(selected[0][i] == full[0][RP_SOL_SPH_R + i - 1U]);
            } else {
                assert(selected[0][i] == full[0][RP_SOL_SPH_VR + i - 4U]);
            }
        }
    }
    assert(rp_solution_reconstruct_canonical_family_batch(0.5,
        &states[0][0], 3U, RP_SOLUTION_FAMILY_ELEMENTS, &elements[0][0],
        statuses) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    for (i = 0U; i < 6U; ++i) {
        assert(elements[0][i] == full[0][RP_SOL_SEMIMAJOR + i]);
        assert(elements[1][i] == 0.0);
    }
    assert(memcmp(original, states, sizeof(states)) == 0);
    assert(rp_solution_reconstruct_canonical_family_batch(NAN, NULL, 0U,
        RP_SOLUTION_FAMILY_ELEMENTS, NULL, NULL) == RP_KERR_STATUS_OK);
    assert(rp_solution_reconstruct_canonical_family_batch(NAN, NULL, 1U,
        RP_SOLUTION_FAMILY_ELEMENTS, &elements[0][0], statuses)
        == RP_KERR_STATUS_NULL_POINTER);
    assert(statuses[0] == RP_KERR_STATUS_NULL_POINTER);
    assert(rp_solution_reconstruct_canonical_family_batch(0.5,
        &states[0][0], 1U, (rp_solution_family)99, &selected[0][0],
        statuses) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(statuses[0] == RP_KERR_STATUS_INVALID_PARAMETER);
    for (i = 0U; i < 7U; ++i) assert(selected[0][i] == 0.0);
    assert(rp_solution_reconstruct_canonical_family_batch(2.0,
        &states[0][0], 1U, RP_SOLUTION_FAMILY_ELEMENTS, &elements[0][0],
        statuses) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(statuses[0] == RP_KERR_STATUS_INVALID_PARAMETER);
    return 0;
}
