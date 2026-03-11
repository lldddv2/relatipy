from astropy import units as u
from ..constants import _c_SI, _G_SI, _M_ref_kg, _L_ref, _T_ref
from numpy import asarray


class Validator:
    def __init__(self):
        self.associated_units_validation = {
            u.kg: self.validate_mass,
            u.s: self.validate_time,
            u.m: self.validate_length,
            u.rad: self.validate_angle,
            u.m / u.s: self.validate_velocity,
            u.rad / u.s: self.validate_angular_velocity,
        }

    def validate_mass(self, mass):
        """kg → unidades geométricas (M/M_ref, donde M_ref = 1 masa solar)"""
        if isinstance(mass, u.Quantity):
            try:
                return mass.to(u.kg).value / _M_ref_kg
            except:
                raise ValueError(f"Invalid mass units: {mass.unit}")
        if mass < 0:
            raise ValueError("Mass must be positive")
        return mass

    def validate_length(self, length):
        """metros → unidades geométricas (r / L_ref, donde L_ref = GM_ref/c²)"""
        if isinstance(length, u.Quantity):
            try:
                return length.to(u.m).value / _L_ref
            except:
                raise ValueError(f"Invalid length units: {length.unit}")
        return length

    def validate_time(self, time):
        """segundos → unidades geométricas (t / T_ref, donde T_ref = GM_ref/c³)"""
        if isinstance(time, u.Quantity):
            try:
                return time.to(u.s).value / _T_ref
            except:
                raise ValueError(f"Invalid time units: {time.unit}")
        return time

    def validate_angle(self, angle):
        """radianes → adimensional (sin cambio)"""
        if isinstance(angle, u.Quantity):
            try:
                return angle.to(u.rad).value
            except:
                raise ValueError(f"Invalid angle units: {angle.unit}")
        return angle

    def validate_velocity(self, velocity):
        """m/s → v/c (adimensional)"""
        if isinstance(velocity, u.Quantity):
            try:
                return velocity.to(u.m / u.s).value / _c_SI
            except:
                raise ValueError(f"Invalid velocity units: {velocity.unit}")
        if velocity > 1:
            raise ValueError("Velocity must be less than the speed of light")
        return velocity

    def validate_angular_velocity(self, angular_velocity):
        """rad/s → adimensional (r_geom * dθ/dt_geom = r[m]/L_ref * dθ/dt[rad/s] * T_ref)
        Como vs = r_geom * dθ/dt_geom y en G=c=1 L_ref = T_ref (ambos GM/c² y GM/c³ difieren por c=1),
        el factor es simplemente T_ref (= L_ref en unidades c=1).
        """
        if isinstance(angular_velocity, u.Quantity):
            try:
                return angular_velocity.to(u.rad / u.s).value * _T_ref
            except:
                raise ValueError(f"Invalid angular velocity units: {angular_velocity.unit}")
        return angular_velocity

    def validate_scalar(self, scalar):
        if isinstance(scalar, u.Quantity):
            try:
                return self.associated_units_validation[scalar.unit](scalar)
            except:
                raise ValueError(f"Invalid scalar units: {scalar.unit}")
        return scalar

    def validate_vector(self, vector):
        vector_si = [x.si if isinstance(x, u.Quantity) else x for x in vector]
        new_vector = []
        for i in vector_si:
            new_vector.append(self.validate_scalar(i))
        return asarray(new_vector)


validator = Validator()