/* Complete elliptic integrals and Carlson forms against mpmath references. */
#include "force/models/elliptic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_relative(const char *label, double got, double want, double tol)
{
    const double err = isinf(want) ? (got == want ? 0.0 : INFINITY)
                                    : fabs(got - want) / fabs(want);
    if (!(err <= tol)) {
        printf("FAIL %s: got %.17g want %.17g (rel %.2e)\n", label, got, want, err);
        ++failures;
    }
}

/* m (exact double), K(m), E(m), K(m) - E(m): mpmath 1.3.0, 60 digits. */
static const double ke_table[][4] = {
    {-10.0, 0.7908718902387384752, 3.6391380384177681635, -2.8482661481790296883},
    {-1.0, 1.3110287771460599052, 1.910098894513856009, -0.59907011736779610372},
    {0.0, 1.5707963267948966192, 1.5707963267948966192, 0.0},
    {1e-12, 1.5707963267952893183, 1.5707963267945039201, 7.8539816339774281813e-13},
    {1e-06, 1.5707967194941992113, 1.5707959340957412894, 7.8539845792194362595e-7},
    {0.01, 1.5747455615173559528, 1.5668619420216682911, 0.0078836194956876616133},
    {0.1, 1.6124413487202194007, 1.5307576368977632002, 0.081683711822456200482},
    {0.3, 1.7138894481787910555, 1.445363064412665267, 0.26852638376612578857},
    {0.5, 1.8540746773013719184, 1.3506438810476755025, 0.50343079625369641591},
    {0.5000001, 1.8540747620226896236, 1.3506438307045941847, 0.50343093131809543888},
    {0.7, 2.0753631352924690784, 1.2416705679458227773, 0.83369256734664630108},
    {0.9, 2.5780921133481732927, 1.1047747327040733079, 1.4733173806440999848},
    {0.99, 3.6956373629898742386, 1.0159935450252239477, 2.679643817964650291},
    {0.999999, 8.2940514636010622019, 1.000003897026172166, 7.2940475665748900359},
};

/* m1 = 1 - m (exact double), K, E: mpmath, 700 digits. */
static const double ke_complement_table[][3] = {
    {1e-300, 346.77405831022674321, 1.0},
    {1e-30, 35.925070756030575837, 1.0},
    {1e-16, 19.806975105072256572, 1.0000000000000009653},
    {1e-10, 12.899219826387599516, 1.0000000006199609913},
    {2.5e-05, 6.6846472569762091551, 1.0000773083030057037},
};

/* x, y, z, R_F, R_D. */
static const double carlson_table[][5] = {
    {2.0, 3.0, 4.0, 0.58408284167715170669, 0.16510527294261053349},
    {0.0, 0.5, 1.0, 1.8540746773013719184, 3.0205847775221784955},
    {1e-08, 1.0, 3.0, 1.1713623527603245535, 0.66485332456866569351},
    {0.25, 4.0, 7.0, 0.59882474366736952402, 0.12319686458465518183},
};

int main(void)
{
    const double tol = 6.7e-16; /* 3 ulp */
    double k, e, kme;

    for (size_t i = 0; i < sizeof ke_table / sizeof ke_table[0]; ++i) {
        if (rp_elliptic_ke(ke_table[i][0], &k, &e) != RP_ELLIPTIC_STATUS_OK) {
            printf("FAIL status at m=%g\n", ke_table[i][0]);
            ++failures;
            continue;
        }
        check_relative("K(m)", k, ke_table[i][1], tol);
        check_relative("E(m)", e, ke_table[i][2], tol);
        const double m = ke_table[i][0];
        if (m >= 0.0 && m < 1.0) {
            /* The hot-path AGM, valid on [0, 1). */
            rp_elliptic_agm(m, 1.0 - m, &k, &kme);
            check_relative("agm K", k, ke_table[i][1], tol);
            check_relative("agm E", k - kme, ke_table[i][2],
                           m > 0.5 ? 2.0e-15 : tol);
            if (m > 0.0) {
                check_relative("agm K-E", kme, ke_table[i][3], tol);
            }
        }
    }
    for (size_t i = 0; i < sizeof ke_complement_table / sizeof ke_complement_table[0]; ++i) {
        rp_elliptic_ke_complementary(ke_complement_table[i][0], &k, &e);
        check_relative("K(1-m1)", k, ke_complement_table[i][1], tol);
        check_relative("E(1-m1)", e, ke_complement_table[i][2], tol);
        /* The AGM's E = K - (K - E) cancels as m -> 1: about K ulp. */
        rp_elliptic_agm(1.0 - ke_complement_table[i][0], ke_complement_table[i][0], &k, &kme);
        check_relative("agm K(1-m1)", k, ke_complement_table[i][1], tol);
        check_relative("agm E(1-m1)", k - kme, ke_complement_table[i][2],
                       1.0e-15 * ke_complement_table[i][1]);
    }
    /* m = 1 and the small-m relation (K - E) / m -> pi / 4. */
    rp_elliptic_ke(1.0, &k, &e);
    if (!(isinf(k) && e == 1.0)) {
        printf("FAIL m = 1: K=%g E=%g\n", k, e);
        ++failures;
    }
    rp_elliptic_agm(1.0e-20, 1.0, &k, &kme);
    check_relative("(K-E)/m at m=1e-20", kme / 1.0e-20, 0.78539816339744830962, tol);

    for (size_t i = 0; i < sizeof carlson_table / sizeof carlson_table[0]; ++i) {
        const double *row = carlson_table[i];
        check_relative("R_F", rp_carlson_rf(row[0], row[1], row[2]), row[3], tol);
        check_relative("R_F perm", rp_carlson_rf(row[2], row[0], row[1]), row[3], tol);
        check_relative("R_D", rp_carlson_rd(row[0], row[1], row[2]), row[4], tol);
        check_relative("R_D swap", rp_carlson_rd(row[1], row[0], row[2]), row[4], tol);
    }
    check_relative("R_D(0,2,1)", rp_carlson_rd(0.0, 2.0, 1.0), 1.7972103521033883112, tol);

    /* Invalid input. */
    if (rp_elliptic_ke(1.5, &k, &e) != RP_ELLIPTIC_STATUS_INVALID_ARGUMENT
        || rp_elliptic_ke(NAN, &k, &e) != RP_ELLIPTIC_STATUS_INVALID_ARGUMENT
        || rp_elliptic_ke(0.5, NULL, NULL) != RP_ELLIPTIC_STATUS_NULL_POINTER
        || !isnan(rp_carlson_rf(-1.0, 1.0, 1.0))
        || !isnan(rp_carlson_rd(1.0, 1.0, 0.0))
        || !isinf(rp_carlson_rf(0.0, 0.0, 1.0))) {
        printf("FAIL invalid-input handling\n");
        ++failures;
    }

    if (failures != 0) {
        printf("test_thin_disk_elliptic: %d failure(s)\n", failures);
        return EXIT_FAILURE;
    }
    printf("test_thin_disk_elliptic: ok\n");
    return EXIT_SUCCESS;
}
