# OBI-Warp

Ordered Bijective Interpolated Warping (OBI-Warp) aligns matrices along a
single axis using Dynamic Time Warping (DTW) and a one-to-one (bijective)
interpolated warp function.  OBI-Warp harnesses the non-linear, comprehensive
alignment power of DTW and builds on the discrete, non-bijective output of DTW
to give natural interpolants that can be used across multiple datasets.

OBI-Warp was developed specifically for the chromatographic alignment of
complex mass spectrometry (MS) proteomics data.  Using high confidence MS/MS
identifications as time standards, OBI-Warp default parameters have been
optimized to give accurate alignments under a variety of real-world conditions
including datasets with little overlapping signal.  Command-line options to
override defaults are available (e.g., gap penalty, local weights and number
of bijective anchors).  Though developed for MS proteomics data, OBI-Warp is
suited to a wide variety of alignment problems.

Pearson's correlation coefficient, covariance, dot product, and Euclidean
distance have been implemented as the available vector similarity functions.
Redundant calculations for correlation coefficient and covariance are cached
in the n x m comparisons to give the algorithmic equivalent of calculating the
dot product.

The dynamic programming algorithm is written to allow any arbitrary gap
penalty function, or users may use a linear initiation and elongation penalty.
Local weighting schemes may also be controlled.

### Links

* [Project Summary Page](://sourceforge.net/projects/obi-warp/)
* [Download OBI-Warp](://sourceforge.net/project/showfiles.php?group_id=161548)
* [Analytical Chemistry Publication](://dx.doi.org/10.1021/ac0605344)
* [Supplementary Material for Publication](://bioinformatics.icmb.utexas.edu/obi-warp/) (*this server is currently down*)

## Building

Building is tested rigorously on Ubuntu and should work fine on any POSIX
system.  Windows compilation under VC++ should be possible and compilation
with cygwin, or msys/mingw should work without any problems.

### Prerequisites

[ruby](http://www.ruby-lang.org) and *rake*. *rake* comes standard with the
newer ruby (1.9.X), but it can also be installed using
[rubygems](http://rubygems.org/pages/download) (for ruby versions below 1.9):

    gem install rake

Optional: *valgrind* can be used for memory testing and
[matrix2png](http://www.bioinformatics.ubc.ca/matrix2png/download.html) can be
used to create images of the various matrices.

### Compiling

    rake
    # creates bin/obiwarp   (notice it is in the *bin* directory)

This will compile all the code and link it into the obiwarp executable.
**NOTE**: All executables will be created in the **bin** directory (testing
executables stay in the lib dir).

If you want to explore other options:

    rake -T

From within the *lib* directory, the generation of any file, intermediate, or
test can be invoked by name.  For example:

    # build the lmat2lmata binary:
    rake lmat2lmata    # note, it will be in created in the bin directory
    # this also works
    rake ../bin/lmat2lmata

From the top directory, all the tasks (not necessarily filetasks) are available:

    # from the top level directory
    rake memtest
    # from the lib directory
    rake memtest

### Installation

Binaries are compiled and depositied in the **bin** folder. System
installation is left to the user, but it can be as simple as:

    rake      # make sure it is compiled and linked
    sudo cp bin/obiwarp /usr/local/bin/

### Testing

    # run all tests
    rake test
    # if you have valgrind installed:
    rake memtest
    rake test_cmdparser.rb  # for tests written in ruby
    rake run_test_dynprog   # for tests written in cxxtest 'run_test_<whatever>'


tasks found in the lib Rakefile are also available from the top level Rakefile using
the 'lib:' prefix:

    rake clean      # from inside the lib directory
    rake lib:clean  # from the top dir

