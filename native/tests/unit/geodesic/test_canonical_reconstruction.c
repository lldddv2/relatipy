#include "geodesic/solution/reconstruct.h"
#include "geodesic/initial/convert.h"

#include <assert.h>
#include <math.h>
#include <string.h>

static void preserves_integrated_velocity(void)
{
    const double bl[7] = {1.0, 8.0, 1.1, 20.0, -0.01, 0.002, 0.02};
    double state[8];
    double before[8];
    double rows[29];
    double cartesian[7];
    double normalized[29];
    rp_kerr_status status;
    size_t i;
    assert(rp_initial_from_bl(0.5, bl, state) == RP_KERR_STATUS_OK);
    /* Deliberate normalization drift must be visible after conversion. */
    for (i = 4U; i < 8U; ++i) {
        state[i] *= 1.01;
    }
    memcpy(before, state, sizeof(state));
    assert(rp_solution_reconstruct_canonical_batch(
        0.5, state, 1U, rows, cartesian, &status
    ) == RP_KERR_STATUS_OK);
    assert(status == RP_KERR_STATUS_OK);
    assert(memcmp(state, before, sizeof(state)) == 0);
    for (i = 0U; i < 4U; ++i) {
        assert(rows[i] == state[i]);
        assert(rows[i + 7U] == state[i + 4U]);
    }
    assert(rows[RP_SOL_SPH_PHI] == state[3]);
    for (i = 0U; i < 3U; ++i) {
        assert(rows[RP_SOL_BL_VR + i] == state[5U + i] / state[4]);
        assert(rows[RP_SOL_UX + i] == state[4] * cartesian[4U + i]);
        assert(rows[RP_SOL_SPH_UR + i] == state[4] * rows[RP_SOL_SPH_VR + i]);
    }
    /* Cartesian interpolation retains its independent normalization path. */
    assert(rp_solution_reconstruct_batch(
        0.5, cartesian, 1U, normalized, &status
    ) == RP_KERR_STATUS_OK);
    assert(fabs(rows[RP_SOL_UT] / normalized[RP_SOL_UT] - 1.01) < 1e-13);
}

static void preserves_high_condition_partial_state(void)
{
    const double state[8] = {0.1, 2.0 + 1e-8, 1.1, 0.2, 2e8, -1.0, 0.0, 0.0};
    double rows[29];
    double cartesian[7];
    rp_kerr_status status;
    size_t i;
    assert(rp_solution_reconstruct_canonical_batch(
        0.0, state, 1U, rows, cartesian, &status
    ) == RP_KERR_STATUS_OK);
    for (i = 0U; i < RP_SOL_SEMIMAJOR; ++i) {
        assert(isfinite(rows[i]));
    }
    for (i = 0U; i < 4U; ++i) {
        assert(rows[i] == state[i]);
        assert(rows[i + 7U] == state[i + 4U]);
    }
}

static void clears_failed_rows_and_continues(void)
{
    double states[3][8] = {
        {0, 8, 1.1, 0, 1.2, 0, 0, .02},
        {0, 8, 0, 0, 1.2, 0, 0, .02},
        {0, 8, 1.1, 0, 1.2, 0, 0, .02}
    };
    double rows[3][29];
    double cartesian[3][7];
    rp_kerr_status statuses[3];
    size_t i;
    assert(rp_solution_reconstruct_canonical_batch(
        0, &states[0][0], 3U, &rows[0][0], &cartesian[0][0], statuses
    ) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert(statuses[0] == RP_KERR_STATUS_OK);
    assert(statuses[1] == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert(statuses[2] == RP_KERR_STATUS_OK);
    for (i = 0U; i < 29U; ++i) assert(rows[1][i] == 0.0);
    for (i = 0U; i < 7U; ++i) assert(cartesian[1][i] == 0.0);
    states[0][4] = -1;
    assert(rp_solution_reconstruct_canonical_batch(
        0, &states[0][0], 1U, &rows[0][0], &cartesian[0][0], statuses
    ) == RP_KERR_STATUS_NON_TIMELIKE_VELOCITY);
    states[0][4] = NAN;
    assert(rp_solution_reconstruct_canonical_batch(
        0, &states[0][0], 1U, &rows[0][0], &cartesian[0][0], statuses
    ) == RP_KERR_STATUS_NONFINITE_INPUT);
    assert(rp_solution_reconstruct_canonical_batch(
        0, NULL, 1U, &rows[0][0], &cartesian[0][0], statuses
    ) == RP_KERR_STATUS_NULL_POINTER);
    for (i = 0U; i < 29U; ++i) assert(rows[0][i] == 0.0);
    assert(rp_solution_reconstruct_canonical_batch(
        NAN, NULL, 0U, NULL, NULL, NULL
    ) == RP_KERR_STATUS_OK);
}

int main(void)
{
    preserves_integrated_velocity();
    preserves_high_condition_partial_state();
    clears_failed_rows_and_continues();
    return 0;
}
