# DAD-drift

The Droplet and Atmospheric Dispersion (DAD) drift model, is new modeling approach that combines a mechanistic droplet model, a micrometeorological model, and a three-dimensional Gaussian diffusion model. It was developed for spray drift prediction after ground application using boom sprayers and allows for drift curve prediction as well as spray drift deposition prediction at the landscape scale.

## License

DAD drift is released under GPL-2.0-only, see [license text](LICENSE.md) for more details.

## Description

This repository contains the [theoretical documentation](documentation/), some [example projects](examples/), the [source code](source_code/), and the model [executable](executable/) of the `DAD drift` model.

### Theoreticcal Documentation

The [theoretical documentation](documentation/) describes the theoretical background as well as the input/output format of the `DAD drift` model in great detail. The purpose of this documentation is to clearly communicate all physical concepts, empirical equations and modeling principals at the core of `DAD drift` and found the model within the scientific literature by citing sources for all equation when sensible. This documentation should help a potential model user to better understand formats and units of model inputs and enable them to interpret model outputs. Furthermore, we want to open the underlying mechanics of `DAD drift` up to critical review of our peers.

### Example Projects

The [example projects](examples/) contain four exemplary use cases of `DAD drift`, with all necessary input files populated and the correct folder structure. After downloading and unpacking the projects a simple double click on the executable will run the project.
1. A [project](examples/drift_curve) where the model is used for `drift curve prediction`.
2. A [project](examples/single_field) where the model is used for `landscape-level drift prediction` for a single field.
3. A [project](examples/multi_field_1) where the model is used for `landscape-level drift prediction` for multiple fields, where all fields are contained in one input file.
4. A [project](examples/multi_field_2) where the model is used for `landscape-level drift prediction` for multiple fields, where each field is contained in a separate input file.

### Source Code

`DAD-drift` is written in `Fortran 90`, the [source code](source_code/) contains all `source files` of the model as well as a `Makefile` for compilation of the executable.

### Executable

A pre-compiled model [executable](executable/) is provided for download.
