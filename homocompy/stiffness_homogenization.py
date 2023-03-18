import numpy as np


def rule_of_mixtures(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v):
    """
    Calculates effective composite properties using the rule of mixtures:
    https://en.wikipedia.org/wiki/Rule_of_mixtures
        -> upper: E = fvf * Ef + (1-fvf) * Em
        -> lower: E = fvf * Em / (Vf * Em + (1-fvf) * Ef)

    The fibre material is assumed to be transversally isotropic.
    The matrix material is assumed to be isotropic

    :param fvf: float or array_like representing the fibre volume fractions to evaluate
    :param f_e11: float, young's modulus 11 of fibre material
    :param f_e22: float, young's modulus 22 of fibre material
    :param f_v12: float, poisson coefficient 12 of fibre material
    :param f_g12: float, shear modulus 12 of fibre material
    :param f_g23: float, shear modulus 23 of fibre material
    :param m_e: float, young's modulus of matrix material, assumed isotropic
    :param m_v: float, poisson coefficient of matrix material, assumed isotropic
    :return: A tuple of effective composite properties e11, e22, v12, v23, g12 and g23
    """

    mvf = 1 - fvf  # matrix volume fraction
    m_g = _calc_isotropic_shear_mod(m_e, m_v)

    e11 = f_e11 * fvf + mvf * m_e
    v12 = f_v12 * fvf + mvf * m_v

    e22 = m_e * f_e22 / (m_e * fvf + mvf * f_e22)
    g12 = m_g * f_g12 / (m_g * fvf + mvf * f_g12)
    g23 = m_g * f_g23 / (m_g * fvf + mvf * f_g23)

    v23 = e22 / (2 * g23) - 1

    return e11, e22, v12, v23, g12, g23


def chamis_model(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v):
    """
    Calculates effective composite properties using Chamis' model.

    The fibre material is assumed to be transversally isotropic.
    The matrix material is assumed to be isotropic

    :param fvf: float or array_like representing the fibre volume fractions to evaluate
    :param f_e11: float, young's modulus 11 of fibre material
    :param f_e22: float, young's modulus 22 of fibre material
    :param f_v12: float, poisson coefficient 12 of fibre material
    :param f_g12: float, shear modulus 12 of fibre material
    :param f_g23: float, shear modulus 23 of fibre material
    :param m_e: float, young's modulus of matrix material, assumed isotropic
    :param m_v: float, poisson coefficient of matrix material, assumed isotropic
    :return: A tuple of effective composite properties e11, e22, v12, v23, g12 and g23

    References
    ----------
    .. [1] C. C. Chamis. Mechanics of Composite Materials: Past, Present and Future. Technical report, NASA Lewis
           Research Center Cleveland, Ohio, Blacksburg, Virginia, 1984.
    """
    m_g = _calc_isotropic_shear_mod(m_e, m_v)

    e11 = f_e11 * fvf + (1 - fvf) * m_e
    v12 = f_v12 * fvf + (1 - fvf) * m_v

    e22 = m_e / (1 - np.sqrt(fvf) * (1 - m_e / f_e22))
    g12 = m_g / (1 - np.sqrt(fvf) * (1 - m_g / f_g12))
    g23 = m_g / (1 - np.sqrt(fvf) * (1 - m_g / f_g23))

    v23 = e22 / (2 * g23) - 1

    return e11, e22, v12, v23, g12, g23


def mori_tanaka(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v):
    """
    Calculates effective composite properties using a closed-form expression of the mori-tanaka model.

    The fibre material is assumed to be transversally isotropic.
    The matrix material is assumed to be isotropic

    :param fvf: float or array_like representing the fibre volume fractions to evaluate
    :param f_e11: float, young's modulus 11 of fibre material
    :param f_e22: float, young's modulus 22 of fibre material
    :param f_v12: float, poisson coefficient 12 of fibre material
    :param f_g12: float, shear modulus 12 of fibre material
    :param f_g23: float, shear modulus 23 of fibre material
    :param m_e: float, young's modulus of matrix material, assumed isotropic
    :param m_v: float, poisson coefficient of matrix material, assumed isotropic
    :return: A tuple of effective composite properties e11, e22, v12, v23, g12 and g23

    References
    ----------
    .. [1] S. G. Abaimov, A. A. Khudyakova, S. V. Lomov, On the closed form expression of the mori–tanaka theory
           prediction for the engineering constants of a unidirectional fiber-reinforced ply, Composite structures
           142 (2016) 1–6
    """

    f_v23 = f_e22 / (2 * f_g23) - 1  # calculate 23 poisson's coefficient of fibre material

    z1 = np.power(
        -2 * (1 - fvf) * np.square(f_v12) / f_e11
        + (1 - fvf) * (1 - f_v23) / f_e22
        + ((1 + m_v) * (1 + fvf * (1 - 2 * m_v))) / m_e,
        -1,
    )
    z2 = f_e22 * (3 + fvf - 4 * m_v) * (1 + m_v) + (1 - fvf) * m_e * (1 + f_v23)

    e11 = (fvf * f_e11 + (1 - fvf) * m_e + 2 * fvf * (1 - fvf) * z1 * np.square(f_v12 - m_v))

    e22 = (e11 / (1 - np.square(m_v))) / (
            1 / (1 - np.square(m_v))
            + 2 * fvf * (e11 / z2) * (1 + f_v23 - f_e22 / m_e * (1 + m_v))
            + fvf
            * z1
            * (f_e11 / m_e)
            * ((1 + m_v) / m_e - 2 / f_e11 + (1 - f_v23) / f_e22)
    )

    v12 = m_v + 2 * fvf * z1 / m_e * (f_v12 - m_v) * (1 - np.square(m_v))

    g12 = (
        m_e
        / (2 * (1 - fvf) * (1 + m_v))
        * (
                1
                + fvf
                - (4 * fvf)
                / (1 + fvf + 2 * (1 - fvf) * (f_g12 / m_e) * (1 + m_v))
        )
    )
    g23 = m_e * np.power(
        2 * (1 + m_v)
        + fvf
        / (
                (1 - fvf) / (8 * (1 - np.square(m_v)))
                + f_g23 / (m_e - 2 * f_g23 * (1 + m_v))
        ),
        -1,
    )

    v23 = e22 / (2 * g23) - 1

    return e11, e22, v12, v23, g12, g23


def _calc_isotropic_shear_mod(e: float, v: float) -> float:
    """
    Calculates the shear modulus of an isotropic material based on its young's modulus and poisson coefficient.

    :param e: float, young's modulus of the material
    :param v: float, poisson coefficient of the material
    :return: float, shear modulus of the material
    """
    return e / (2 * (v + 1))