import numpy as np
import matplotlib.pyplot as plt
from HomoComPy.stiffness_homogenization import rule_of_mixtures, chamis_model, mori_tanaka


def example():
    """
    This example plots the homogenized effective composite properties over a range of fibre volume fractions,
    using several homogenization models.
    """

    # Fibre properties - units in GPa - transversally isotropic
    f_e11 = 230.0  # Young's modulus 11
    f_e22 = 20.0  # Young's modulus 22
    f_v12 = 0.2  # Poisson coefficient 12
    f_g12 = 30.0  # Shear modulus 12
    f_g23 = 7.0  # Shear modulus 23

    # Matrix properties - units in GPa - isotropic
    m_e = 4.0  # Young's modulus
    m_v = 0.3  # Poisson coefficient

    # Fibre fraction array: 0% to 100% in steps of 1%
    fvf = np.linspace(0.01, 0.99, 100)

    # Calculate effective properties
    # Results is a tuple of floats: (e11, e22/e33, v12/v13, v23, g12/g13, g23)
    results_rom = rule_of_mixtures(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v)
    results_chamis = chamis_model(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v)
    results_mori_tanaka = mori_tanaka(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v)

    plt_id = 1
    prop_name = ["E11 [GPa]", "E22/E33 [GPa]", "v12/v13 [-]", "v23 [-]", "G12/G13 [GPa]", "G23 [GPa]"][plt_id]

    plt.plot(fvf * 100, results_rom[plt_id], label="ROM")
    plt.plot(fvf * 100, results_chamis[plt_id], label="Chamis")
    plt.plot(fvf * 100, results_mori_tanaka[plt_id], label="Mori-Tanaka")
    plt.legend()
    plt.xlabel("Fibre Volume Fraction [%]")
    plt.ylabel(f"{prop_name}")
    plt.title("Homogenized Effective Property")
    plt.show()

if __name__ == "__main__":
    example()