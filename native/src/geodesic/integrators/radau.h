#ifndef RELATIPY_GEODESIC_INTEGRATORS_RADAU_H
#define RELATIPY_GEODESIC_INTEGRATORS_RADAU_H

#include "integrator.h"

rp_integrator_status rp_radau_integrate(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double initial_independent_variable,
    double final_independent_variable,
    double state[],
    rp_integrator_stats *stats
);

#endif /* RELATIPY_GEODESIC_INTEGRATORS_RADAU_H */
