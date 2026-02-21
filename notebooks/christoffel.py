from numpy import cos, sin, tan
import numpy

def get_christoffel_symbols(xs, a, R_s):
        """
        Returns the Christoffel symbols of the metric.

        Parameters
        ----------
        xs : list
            List of coordinates [x0, x1, x2, x3].
        """
        x0, x1, x2, x3 = xs

        ############### Auxiliary variables ##############

        cos_I2_x2I_ = cos(2 * x2)
        cos_I4_x2I_ = cos(4 * x2)
        cos_Ix2I_ = cos(x2)
        sin_I2_x2I_ = sin(2 * x2)
        sin_Ix2I_ = sin(x2)
        tan_Ix2I_ = tan(x2)

        R_s_pow2 = R_s**2
        _1_over_tan_Ix2I_ = 1 / tan_Ix2I_
        a_pow2 = a**2
        a_pow4 = a**4
        a_pow6 = a**6
        cos_Ix2I__pow2 = cos_Ix2I_**2
        cos_Ix2I__pow4 = cos_Ix2I_**4
        sin_Ix2I__pow2 = sin_Ix2I_**2
        sin_Ix2I__pow4 = sin_Ix2I_**4
        x1_pow2 = x1**2
        x1_pow3 = x1**3
        x1_pow4 = x1**4
        x1_pow5 = x1**5
        x1_pow6 = x1**6

        # Modified manually for optimization
        a_pow2_cos_Ix2I__pow2 = a_pow2 * cos_Ix2I__pow2
        R_s_a_pow2_cos_Ix2I__pow2_x1 = R_s * a_pow2_cos_Ix2I__pow2 * x1
        R_s_x1_pow3 = R_s * x1_pow3
        _2_x1_pow2 = 2 * x1_pow2
        _R_s_x1 = -R_s * x1
        _a_pow2_cos_Ix2I__pow2_x1_pow2 = -a_pow2_cos_Ix2I__pow2 * x1_pow2
        a_pow2_cos_I2_x2I_ = a_pow2 * cos_I2_x2I_
        _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_ = (
            -a_pow2 * sin_I2_x2I_ / (a_pow2_cos_I2_x2I_ + a_pow2 + _2_x1_pow2)
        )
        a_pow2_x1_pow2 = a_pow2 * x1_pow2
        _a_pow2_x1_pow2 = -a_pow2_x1_pow2
        _cos_I4_x2I_ = -cos_I4_x2I_
        _x1_pow4 = -x1_pow4
        a_pow4_cos_Ix2I__pow2 = a_pow4 * cos_Ix2I__pow2

        ##################################################

        ############### Christoffel symbols ##############
        Gamma = numpy.zeros((4, 4, 4))

        # Only non-zero components and symmetry properties

        # mu = 0
        Gamma[0, 0, 1] = (
            R_s
            * (
                a_pow2_x1_pow2 * sin_Ix2I__pow2
                + a_pow4 * sin_Ix2I__pow2
                - a_pow4
                + x1_pow4
            )
            / (
                2
                * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
                * (_R_s_x1 + a_pow2 + x1_pow2)
            )
        )
        Gamma[0, 0, 2] = (
            2
            * R_s
            * _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_
            * x1
            / (_2_x1_pow2 + a_pow2 + a_pow2_cos_I2_x2I_)
        )
        Gamma[0, 1, 3] = (
            R_s
            * a
            * sin_Ix2I__pow2
            * (
                _a_pow2_cos_Ix2I__pow2_x1_pow2
                + _a_pow2_x1_pow2
                + 3 * _x1_pow4
                + a_pow4_cos_Ix2I__pow2
            )
            / (
                2
                * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
                * (_R_s_x1 + a_pow2 + x1_pow2)
            )
        )
        Gamma[0, 2, 3] = (
            R_s
            * a**3
            * cos_Ix2I_
            * sin_Ix2I_**3
            * x1
            / (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
        )

        Gamma[0, 1, 0] = Gamma[0, 0, 1]
        Gamma[0, 2, 0] = Gamma[0, 0, 2]
        Gamma[0, 3, 1] = Gamma[0, 1, 3]
        Gamma[0, 3, 2] = Gamma[0, 2, 3]

        # mu = 1
        Gamma[1, 0, 0] = (
            R_s
            * (
                R_s_a_pow2_cos_Ix2I__pow2_x1
                - R_s_x1_pow3
                + _a_pow2_cos_Ix2I__pow2_x1_pow2
                + _a_pow2_x1_pow2
                + a_pow4_cos_Ix2I__pow2
                + x1_pow4
            )
            / (2 * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 3)
        )
        Gamma[1, 0, 3] = (
            R_s
            * a
            * sin_Ix2I__pow2
            * (
                -R_s_a_pow2_cos_Ix2I__pow2_x1
                + R_s_x1_pow3
                + _x1_pow4
                + a_pow2_cos_Ix2I__pow2 * x1_pow2
                + a_pow2_x1_pow2
                - a_pow4_cos_Ix2I__pow2
            )
            / (2 * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 3)
        )
        Gamma[1, 1, 1] = (
            -R_s * a_pow2_cos_Ix2I__pow2 / 2
            + R_s * x1_pow2 / 2
            + a_pow2 * x1
            + a_pow2_cos_Ix2I__pow2 * x1
        ) / (
            R_s_a_pow2_cos_Ix2I__pow2_x1
            + R_s_x1_pow3
            + _a_pow2_cos_Ix2I__pow2_x1_pow2
            + _x1_pow4
            + a_pow2_x1_pow2
            + a_pow4_cos_Ix2I__pow2
        )
        Gamma[1, 1, 2] = (
            _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_
        )
        Gamma[1, 2, 2] = (
            x1 * (R_s * x1 + a_pow2 - x1_pow2) / (a_pow2_cos_Ix2I__pow2 + x1_pow2)
        )
        Gamma[1, 3, 3] = (
            sin_Ix2I__pow2
            * (
                R_s * a_pow2 * sin_Ix2I__pow2 * x1_pow4 / 2
                + 2 * R_s * a_pow2_cos_Ix2I__pow2 * x1_pow4
                + R_s * a_pow4 * cos_Ix2I__pow4 * x1_pow2
                - R_s * a_pow4 * sin_Ix2I__pow2 * x1_pow2 / 2
                - R_s * a_pow4 * x1_pow2 * (_cos_I4_x2I_ + 1) / 16
                + R_s * a_pow6 * (_cos_I4_x2I_ + 1) / 16
                + R_s * x1_pow6
                - R_s_pow2 * a_pow2 * sin_Ix2I__pow2 * x1_pow3 / 2
                + R_s_pow2 * a_pow4 * x1 * (_cos_I4_x2I_ + 1) / 16
                + a_pow2 * x1_pow5
                - 2 * a_pow2_cos_Ix2I__pow2 * x1_pow5
                - a_pow4 * cos_Ix2I__pow4 * x1_pow3
                + 2 * a_pow4_cos_Ix2I__pow2 * x1_pow3
                + a_pow6 * cos_Ix2I__pow4 * x1
                - x1**7
            )
            / (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 3
        )

        Gamma[1, 2, 1] = Gamma[1, 1, 2]
        Gamma[1, 3, 0] = Gamma[1, 0, 3]

        # mu = 2
        Gamma[2, 0, 0] = (
            4
            * R_s
            * _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_
            * x1
            / (_2_x1_pow2 + a_pow2 + a_pow2_cos_I2_x2I_) ** 2
        )
        Gamma[2, 0, 3] = (
            4
            * R_s
            * a
            * sin_I2_x2I_
            * x1
            * (a_pow2 + x1_pow2)
            / (_2_x1_pow2 + a_pow2 + a_pow2_cos_I2_x2I_) ** 3
        )
        Gamma[2, 1, 1] = (
            -a_pow2
            * cos_Ix2I_
            * sin_Ix2I_
            / (
                R_s_a_pow2_cos_Ix2I__pow2_x1
                + R_s_x1_pow3
                + _a_pow2_cos_Ix2I__pow2_x1_pow2
                - _a_pow2_x1_pow2
                + _x1_pow4
                + a_pow4_cos_Ix2I__pow2
            )
        )
        Gamma[2, 1, 2] = x1 / (a_pow2_cos_Ix2I__pow2 + x1_pow2)
        Gamma[2, 2, 2] = (
            _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_
        )
        Gamma[2, 3, 3] = (
            cos_Ix2I_
            * sin_Ix2I_
            * (
                -2 * R_s_x1_pow3 * a_pow2 * sin_Ix2I__pow2
                - _2_x1_pow2 * a_pow4_cos_Ix2I__pow2
                + _R_s_x1 * a_pow4 * sin_Ix2I__pow4
                + _R_s_x1 * a_pow4 * (_cos_I4_x2I_ + 1) / 4
                + _x1_pow4 * a_pow2
                + 2 * _x1_pow4 * a_pow2_cos_Ix2I__pow2
                - a_pow4 * cos_Ix2I__pow4 * x1_pow2
                - a_pow6 * cos_Ix2I__pow4
                - x1_pow6
            )
            / (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 3
        )

        Gamma[2, 2, 1] = Gamma[2, 1, 2]
        Gamma[2, 3, 0] = Gamma[2, 0, 3]

        # mu = 3
        Gamma[3, 0, 1] = (
            R_s
            * a
            * (-a_pow2_cos_Ix2I__pow2 + x1_pow2)
            / (
                2
                * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
                * (_R_s_x1 + a_pow2 + x1_pow2)
            )
        )
        Gamma[3, 0, 2] = (
            _1_over_tan_Ix2I_ * _R_s_x1 * a / (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
        )
        Gamma[3, 1, 3] = (
            R_s * _a_pow2_cos_Ix2I__pow2_x1_pow2
            + R_s * _a_pow2_x1_pow2 * sin_Ix2I__pow2 / 2
            + R_s * _x1_pow4
            + R_s * a_pow4 * (_cos_I4_x2I_ + 1) / 16
            + 2 * a_pow2_cos_Ix2I__pow2 * x1_pow3
            + a_pow4 * cos_Ix2I__pow4 * x1
            + x1_pow5
        ) / ((a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2 * (_R_s_x1 + a_pow2 + x1_pow2))
        Gamma[3, 2, 3] = (
            _1_over_tan_Ix2I_
            * (
                R_s * a_pow2 * sin_Ix2I__pow2 * x1
                + _2_x1_pow2 * a_pow2
                + 2 * _a_pow2_x1_pow2 * sin_Ix2I__pow2
                - 2 * a_pow4 * sin_Ix2I__pow2
                + a_pow4 * sin_Ix2I__pow4
                + a_pow4
                + x1_pow4
            )
            / (_2_x1_pow2 * a_pow2_cos_Ix2I__pow2 + a_pow4 * cos_Ix2I__pow4 + x1_pow4)
        )

        Gamma[3, 1, 0] = Gamma[3, 0, 1]
        Gamma[3, 2, 0] = Gamma[3, 0, 2]
        Gamma[3, 3, 1] = Gamma[3, 1, 3]
        Gamma[3, 3, 2] = Gamma[3, 2, 3]

        ##################################################

        return Gamma
