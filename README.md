# HomoComPy
This packages helps in computing effective stiffness properties of continuous fibre-matrix composite materials, based on analytical homogenization models. It assumes the fibres are transversally isotropic and the matrix is isotropic

Currently, these models are supported:
- *Rule of Mixtures*
- *Chamis' model*
- *Mori-Tanaka* (based on a closed-form expression)
### Installation
`pip install HomoComPy`

### Get started
How to calculate the effective properties of a fibre / matrix composite with this package:
```
from homocompy import rule_of_mixtures, chamis_model, mori_tanaka

fvf = 0.5  # fibre volume fraction [-]

# Fibre properties - units in GPa - transversally isotropic
f_e11 = 230.0  # Young's modulus 11
f_e22 = 20.0  # Young's modulus 22
f_v12 = 0.2  # Poisson coefficient 12
f_g12 = 30.0  # Shear modulus 12
f_g23 = 7.0  # Shear modulus 23

# Matrix properties - units in GPa - isotropic
m_e = 4.0  # Young's modulus
m_v = 0.3  # Poisson coefficient

# Calculate effective properties
# Results is a tuple of floats: (e11, e22/e33, v12/v13, v23, g12/g13, g23)
results_rom = rule_of_mixtures(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v)
results_chamis = chamis_model(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v)
results_mori_tanaka = mori_tanaka(fvf, f_e11, f_e22, f_v12, f_g12, f_g23, m_e, m_v)

```
