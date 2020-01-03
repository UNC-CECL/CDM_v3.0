# Coastal Dune Model

Coastal Dune Model that includes updated vegetation and wrack dynamics.


## Install

### Prerequisites

To install the Coastal Dune Model from source, you will need to install some
prerequisite packages.
* A C++ compiler
* CMake
* FFTW

If you are using Anaconda, these can be installed through the *conda* command.

    $ conda create -n cdm --file=requirements.txt
    $ conda activate cdm

This creates a new conda environment, *cdm*, with all of the prequisite packages
installed.

### Compile

To build and install cmake, run

    $ mkdir _build && cd _build
    $ cmake ../src -DCMAKE_INSTALL_PREFIX=<path-to-installation>
    $ make -j4
    $ make install

where `<path-to-installation>` is the base directory in which the Coastal Dune
Model will be installed (`/usr/local` is the default).

## Use

With the Coastal Dune Model installed following the above instructions, you
can run it with the *coastal-dune-model* executable that was installed in
`<path-to-installation>/bin`. If this folder is in your *PATH*,

    $ coastal-dune-model params.par

## References

Biel, R. G., Moore, L. J., and Goldstein, E. B. (submitted). Influence of wrack on dune forma-
tion.Journal of Geophysical Research: Biogeosciences.

Biel, R. G., Moore, L. J., Zinnert, J., and Brown, J. (in prep.). Inhibition or facilitation of
dune development on barriers: The influence of back-beach vegetation on barrier island sta-
ble states.Journal of Geophysical Research: Biogeosciences.

Duran, O. and Moore, L. J. (2013). Vegetation controls on the maximum size of coastal dunes.
Proceedings of the National Academy of Sciences, 110(43):17217–17222.

Duran, O., Parteli, E. J., and Herrmann, H. J. (2010). A continuous model for sand dunes:
Review, new developments and application to barchan dunes and barchan dune fields.Earth
Surface Processes and Landforms, 35(13):1591–1600.

Duran Vinent, O. and Moore, L. J. (2015). Barrier island bistability induced by biophysical
interactions.Nature Clim. Change, 5(2):158–162.

Moore, L. J., Vinent, O. D., and Ruggiero, P. (2016). Vegetation control allows autocyclic
formation of multiple dunes on prograding coasts.Geology, 44(7):559–562.

Weng, W. S., Hunt, J. C. R., Carruthers, D. J., Warren, A., Wiggs, G. F. S., Livingstone, I.,
and Castro, I. (1991). Air flow and sand transport over sand-dunes. In Barndorff-Nielsen,
O. E. and Willetts, B. B., editors,Aeolian Grain Transport, pages 1–22, Vienna. Springer
Vienna.
