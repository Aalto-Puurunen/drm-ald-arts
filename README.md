# DRM-ALD-Arts — Diffusion–Reaction Model for Atomic Layer Deposition implemented by Arts

<a href="https://doi.org/10.5281/zenodo.7844551"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7844551.svg" alt="DOI"></a>

## Project description

This is a MATLAB implementation of a diffusion-reaction model for simulating the conformality of atomic layer deposition in high-aspect-ratio structures. 

This MATLAB script implements a diffusion-reaction model developed by Yanguas-Gil and Elam (Chem. Vap. Deposition **18**, 46-52 (2012), DOI: [10.1002/cvde.201106938](https://doi.org/10.1002/cvde.201106938)). The model uses the Damkoler and precursor excess number to determine the temporal evolution of surface coverage  inside a high-aspect-ratio structure as a function of distance from the structure entrance. While the model is flexible in terms of structure geometry, this implementation was created in the context of rectangular trenches. The implementation was written by Dr. Karsten Arts and published by Jänis Järvilehto by request of Prof. Riikka Puurunen.  

## Usage

The simulation is performed by running model.m. 

The simulation parameters can be changed by modifying the variable values assigned after *INPUT VALUES*. An explanation of the variables and their units is given in the comment at the beginning of the script. 

During simulation, the input sticking coefficient, a back-extracted sticking coefficient and the ratio of these values is printed to the command line. The back-extracted value is calculated using the method described in Arts *et al.* (J. Vac. Sci. Technol., A **37**, 030908 (2019), DOI: [10.1116/1.5093620](https://doi.org/10.1116/1.5093620)). 

After simulation, the script writes the resulting saturation and pressure profiles in terms of distance and distance divided by channel height to an Excel file. The file also contains a record of the input parameters. Additionally, the saturation and pressure profiles are plotted as MATLAB figures. 

## Publications

### v1.0.0

K. Arts, M. Utriainen, R. L. Puurunen, W. M. M. Kessels, and H. C. M. Knoops. **Film Conformality and Extracted Recombination Probabilities of O Atoms during Plasma-Assisted Atomic Layer Deposition of SiO2, TiO2, Al2O3, and HfO2**. *The Journal of Physical Chemistry C* 123 (2019) 27030–35. [https://doi.org/10.1021/acs.jpcc.9b08176](https://doi.org/10.1021/acs.jpcc.9b08176).

K. Arts, V. Vandalon, R. L. Puurunen, M. Utriainen, F. Gao, W. M. M. Kessels, and H. C. M. Knoops. **Sticking Probabilities of H2O and Al(CH3)3 during Atomic Layer Deposition of Al2O3 Extracted from Their Impact on Film Conformality**. *Journal of Vacuum Science & Technology A* 37 (2019) 030908. [https://doi.org/10.1116/1.5093620](https://doi.org/10.1116/1.5093620).

K. Arts, S. Deijkers, R. L. Puurunen, W. M. M. Kessels, and H. C. M. Knoops. **Oxygen Recombination Probability Data for Plasma-Assisted Atomic Layer Deposition of SiO2 and TiO2**. *The Journal of Physical Chemistry C* 125 (2021) 8244–52. [https://doi.org/10.1021/acs.jpcc.1c01505](https://doi.org/10.1021/acs.jpcc.1c01505).

## Citing

Please cite as:

K. Arts, **DRM-ALD-Arts — Diffusion–Reaction Model for Atomic Layer Deposition implemented by Arts**, (2023), *Github repository*, [https://github.com/Aalto-Puurunen/drm-ald-arts](https://github.com/Aalto-Puurunen/drm-ald-arts). [https://doi.org/10.5281/zenodo.7844551](https://doi.org/10.5281/zenodo.7844551).

And the related article:

K. Arts, M. Utriainen, R. L. Puurunen, W. M. M. Kessels, and H. C. M. Knoops. **Film Conformality and Extracted Recombination Probabilities of O Atoms during Plasma-Assisted Atomic Layer Deposition of SiO2, TiO2, Al2O3, and HfO2**. *The Journal of Physical Chemistry C* 123 (2019) 27030–35. [https://doi.org/10.1021/acs.jpcc.9b08176](https://doi.org/10.1021/acs.jpcc.9b08176).

The model was originally presented in:

A. Yanguas-Gil, and J. W. Elam. **Self-Limited Reaction-Diffusion in Nanostructured Substrates: Surface Coverage Dynamics and Analytic Approximations to ALD Saturation Times**. *Chemical Vapor Deposition* 18 (2012) 46–52. [https://doi.org/10.1002/cvde.201106938](https://doi.org/10.1002/cvde.201106938).

## Copyright and license

MIT License

Copyright 2023 (c) Karsten Arts, Eindhoven University of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
