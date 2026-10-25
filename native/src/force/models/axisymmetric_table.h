/**
 * @file axisymmetric_table.h
 * @brief Tabulated Newtonian field of an axisymmetric, plane-symmetric mass.
 *
 * Private to the native force models.  A model (such as the thin disk)
 * supplies the exact field of its unit-mass distribution with unit length
 * scale in its own frame: symmetry axis `e_z`, cylindrical radius `d`, height
 * `z`, spherical radius `r = sqrt(d^2 + z^2)` and `mu = |z| / r`.  The table
 * stores the dimensionless, everywhere smooth combinations
 *
 *     A = -r^3 (a_d / d),   C = -r^3 (a_z / z),   P = -r Phi,
 *
 * which are even in `z`, tend to one far away, need no division by `d` or
 * `z` on the axis or in the plane, and give the Cartesian field directly:
 *
 *     a = -(A xi + (C - A) xi_z e_z) / r^3,      Phi = -P / r.
 *
 * Grid: nonuniform nodes in `r` and in `mu`, refined geometrically towards the
 * given edges (radii with a structure width) and towards the plane `mu = 0`.
 * Each quantity is a tensor-product C2 cubic spline (not-a-knot in `r` and at
 * `mu = 1`; zero slope at `mu = 0`, the exact plane symmetry), stored as
 * bicubic Hermite data (value, d/dr, d/dmu, d2/drdmu) per node, with `A` and
 * `C` interleaved so that one cell reads four 64-byte nodes.  Cells are found
 * with a guide table indexed by the leading bits of the double, so no `log`,
 * `atan2` or search runs per evaluation.
 *
 * Below `radial_min`, `A` and `C` are continued as `r^inner_exponent`, the
 * self-similar or hole behaviour declared by the model.  Beyond
 * `radial_match` the model's multipole series is blended in with a C2
 * quintic weight; beyond `radial_far` only the series is used.
 *
 * A built table is immutable and may be shared between threads.
 */

#ifndef RELATIPY_FORCE_MODELS_AXISYMMETRIC_TABLE_H
#define RELATIPY_FORCE_MODELS_AXISYMMETRIC_TABLE_H

#include <stddef.h>
#include <stdint.h>

#define RP_AXISYMMETRIC_MAX_EDGES 4U
/** Highest even multipole order `k` of the far-field series. */
#define RP_AXISYMMETRIC_MAX_FAR_ORDER 24U
#define RP_AXISYMMETRIC_FAR_TERMS (RP_AXISYMMETRIC_MAX_FAR_ORDER / 2U + 1U)

typedef enum {
    RP_AXISYMMETRIC_STATUS_OK = 0,
    RP_AXISYMMETRIC_STATUS_NULL_POINTER,
    RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT,
    RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE,
    /** The model's field callback failed or returned a non-finite value. */
    RP_AXISYMMETRIC_STATUS_FIELD_FAILURE
} rp_axisymmetric_status;

/**
 * Exact field of the unit-mass, unit-length distribution at `(d, z)`, with
 * `d >= 0` and `z >= 0`.  Writes `field[0] = a_d / d`, `field[1] = a_z / z`
 * (both finite limits on the axis and in the plane) and `field[2] = Phi`.
 * Returns zero on success.
 */
typedef int (*rp_axisymmetric_field_fn)(
    double d, double z, void *context, double field[3]
);

/** Grid and asymptotics of one table; every length in the unit length. */
typedef struct {
    double radial_min;     /**< innermost node, > 0 */
    double radial_match;   /**< far-field blend starts, > radial_min */
    double radial_far;     /**< far-field blend ends; outermost node */
    double radial_step;    /**< largest relative spacing h / r, in (0, 0.5] */
    double angular_step;   /**< largest spacing growth h / mu, in (0, 0.5] */
    /** Spacing at an edge or at the plane, in units of its width. */
    double floor_step;
    size_t edge_count;     /**< <= RP_AXISYMMETRIC_MAX_EDGES */
    double edge_radius[RP_AXISYMMETRIC_MAX_EDGES];
    double edge_width[RP_AXISYMMETRIC_MAX_EDGES]; /**< > 0 */
    double plane_width;    /**< width in mu of the plane structure, > 0 */
    /** `A, C` proportional to `r^inner_exponent` below radial_min, > 0. */
    double inner_exponent;
    /** Even order of the far series, <= RP_AXISYMMETRIC_MAX_FAR_ORDER. */
    size_t far_order;
    /**
     * `P = sum_k r^-k G_k(u)` with `u = mu^2` and
     * `G_k(u) = sum_i far_coefficients[k/2][i] u^i`, `i <= k/2`.
     */
    double far_coefficients[RP_AXISYMMETRIC_FAR_TERMS][RP_AXISYMMETRIC_FAR_TERMS];
} rp_axisymmetric_spec;

/** One nonuniform axis with its cell guide; owned by the table. */
typedef struct {
    size_t count;            /**< nodes, >= 4 */
    double *node;            /**< count increasing values */
    double *inverse_step;    /**< count - 1 values 1 / (node[i+1] - node[i]) */
    uint32_t *guide;         /**< guide_count cell indices */
    size_t guide_count;
    uint64_t guide_base;     /**< key of guide[0] */
} rp_axisymmetric_axis;

/** A built table.  Read-only after rp_axisymmetric_table_build(). */
typedef struct rp_axisymmetric_table {
    rp_axisymmetric_axis radius;
    rp_axisymmetric_axis cosine;   /**< mu = |z| / r on [0, 1] */
    /** Node (i, j) at [8 (j n_r + i)]: A, C, A_r, C_r, A_mu, C_mu, A_rmu, C_rmu. */
    double *field;
    /** Node (i, j) at [4 (j n_r + i)]: P, P_r, P_mu, P_rmu. */
    double *potential;
    double radial_min;
    double radial_match;
    double radial_far;
    double inverse_blend_width;
    double inner_exponent;
    size_t far_terms;
    /** Far series of A, C, P: coefficient [k/2][i] of r^-k u^i. */
    double far_a[RP_AXISYMMETRIC_FAR_TERMS][RP_AXISYMMETRIC_FAR_TERMS];
    double far_c[RP_AXISYMMETRIC_FAR_TERMS][RP_AXISYMMETRIC_FAR_TERMS];
    double far_p[RP_AXISYMMETRIC_FAR_TERMS][RP_AXISYMMETRIC_FAR_TERMS];
    size_t bytes;                  /**< heap bytes owned by the table */
} rp_axisymmetric_table;

/**
 * Build a table by calling `field` once per node.  On success `*table`
 * receives a heap-allocated table that the caller releases with
 * rp_axisymmetric_table_free(); on failure `*table` is set to NULL and
 * nothing is leaked.  `context` is passed through and not retained.
 */
rp_axisymmetric_status rp_axisymmetric_table_build(
    const rp_axisymmetric_spec *spec,
    rp_axisymmetric_field_fn field,
    void *context,
    rp_axisymmetric_table **table
);

/** Release a table from rp_axisymmetric_table_build(); NULL is ignored. */
void rp_axisymmetric_table_free(rp_axisymmetric_table *table);

/**
 * Interpolated `A, C, P` of the unit distribution at `r > 0`, `0 <= mu <= 1`,
 * including the inner continuation and the far series.  `out` holds three
 * doubles.  Inputs are not validated.
 */
void rp_axisymmetric_table_fields(
    const rp_axisymmetric_table *table, double r, double mu, double out[3]
);

/**
 * One placed distribution.  `rotation` maps the black-hole frame to the
 * model frame (`xi = rotation x`; its rows are the model axes, the last one
 * the symmetry axis).  Positions are divided by `length_scale`; the field is
 * multiplied by `mass` (geometric units, black-hole mass = 1).  Fill it with
 * rp_axisymmetric_model_init(), or keep `inverse_length_scale` equal to
 * `1 / length_scale` and `rotation` orthonormal when setting fields.
 * The table is borrowed and must outlive every use of the model.
 */
typedef struct {
    const rp_axisymmetric_table *table;
    double mass;
    double length_scale;
    double inverse_length_scale;
    double rotation[3][3];
} rp_axisymmetric_model;

/**
 * Rotation whose last row is the unit `normal` (renormalized when
 * `| |normal| - 1 | <= 1e-10`, rejected otherwise); the first two rows
 * complete a right-handed orthonormal frame.  `rotation` is caller-owned;
 * it is left unchanged on failure.
 */
rp_axisymmetric_status rp_axisymmetric_rotation_from_normal(
    const double normal[3], double rotation[3][3]
);

/**
 * Fill `model` from a table, a finite mass, a positive finite length scale
 * and a unit symmetry axis in the black-hole frame.  On failure `model` is
 * left unchanged.
 */
rp_axisymmetric_status rp_axisymmetric_model_init(
    rp_axisymmetric_model *model,
    const rp_axisymmetric_table *table,
    double mass,
    double length_scale,
    const double normal[3]
);

/**
 * rp_external_acceleration_fn of a placed model (see
 * `force/external_force.h`); `velocity` is unused.  Returns nonzero for a
 * NULL pointer, a non-finite position or the origin.
 */
int rp_axisymmetric_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
);

/** Potential of a placed model at `position`, with the same failures. */
int rp_axisymmetric_potential(
    const double position[3], const void *model, double *potential
);

#endif /* RELATIPY_FORCE_MODELS_AXISYMMETRIC_TABLE_H */
