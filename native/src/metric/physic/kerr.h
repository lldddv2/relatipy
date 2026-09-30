/**
 * @file kerr.h
 * @brief Private boundary for Kerr physical evaluators.
 *
 * This header is internal to the native implementation.  It is not installed
 * and does not define a public API or stable ABI.
 */

#ifndef RELATIPY_NATIVE_METRIC_PHYSIC_KERR_H
#define RELATIPY_NATIVE_METRIC_PHYSIC_KERR_H

#include "relatipy/kerr_geometry.h"

/** Validated coordinates and reusable Kerr point quantities. */
struct rp_kerr_point {
    double x[RP_KERR_DIM];
    double sin_theta;
    double cos_theta;
    double sigma;
    double delta;
};

/** Complete the physical quantities of a point whose inputs were validated. */
rp_kerr_status rp_kerr_physic_prepare_point(
    double mass,
    double spin,
    struct rp_kerr_point *point
);

/** Evaluate the covariant Kerr metric at a prepared point. */
void rp_kerr_physic_evaluate_metric(
    double mass,
    double spin,
    const struct rp_kerr_point *point,
    double metric[RP_KERR_DIM][RP_KERR_DIM]
);

/** Evaluate the contravariant Kerr metric at a prepared point. */
rp_kerr_status rp_kerr_physic_evaluate_inverse_metric(
    double mass,
    double spin,
    const struct rp_kerr_point *point,
    double inverse_metric[RP_KERR_DIM][RP_KERR_DIM]
);

/** Evaluate the frozen-legacy-compatible Kerr connection. */
void rp_kerr_physic_evaluate_christoffel(
    double mass,
    double spin,
    const struct rp_kerr_point *point,
    double christoffel[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM]
);

/** Convert validated coordinate velocities into a four-velocity. */
rp_kerr_status rp_kerr_physic_evaluate_four_velocity(
    double mass,
    double spin,
    const struct rp_kerr_point *point,
    const double coordinate_velocity[3],
    double four_velocity[RP_KERR_DIM]
);

#endif /* RELATIPY_NATIVE_METRIC_PHYSIC_KERR_H */
