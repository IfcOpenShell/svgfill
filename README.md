svgfill
=======

An application to fill areas bounded by unconnected lines in SVG.

Dependencies
------------

* [CGAL 2D Arrangements](https://doc.cgal.org/latest/Arrangement_on_surface_2/index.html) GPL
* [SVG++](http://svgpp.org/) Boost software license

Compilation
-----------

Installation is shown based on the IfcOpenShell windows build script directory output:

    mkdir build
    cd build
    set IFCOPENSHELL_ROOT=..\path\to\ifcopenshell\directory\
    cmake -DBOOST_ROOT=%IFCOPENSHELL_ROOT%\deps\boost_1_67_0                                                    ^
          -DBOOST_LIBRARYDIR=%IFCOPENSHELL_ROOT%\deps\boost_1_67_0\stage\vs2017-Win32\lib                       ^
          -DLIBXML2_INCLUDE_DIR=%IFCOPENSHELL_ROOT%\deps\OpenCOLLADA\Externals\LibXML\include                   ^
          -DLIBXML2_LIBRARIES=%IFCOPENSHELL_ROOT%\deps-vs2017-x86-installed\OpenCOLLADA\lib\opencollada\xml.lib ^
          -DCGAL_INCLUDE_DIR=%IFCOPENSHELL_ROOT%\deps-vs2017-x86-installed\cgal\include                         ^
          -DCGAL_LIBRARY_DIR=%IFCOPENSHELL_ROOT%\deps-vs2017-x86-installed\cgal\lib                             ^
          -DGMP_INCLUDE_DIR=%IFCOPENSHELL_ROOT%\deps-vs2017-x86-installed\mpir                                  ^
          -DGMP_LIBRARY_DIR=%IFCOPENSHELL_ROOT%\deps-vs2017-x86-installed\mpir                                  ^
          -DMPFR_INCLUDE_DIR=%IFCOPENSHELL_ROOT%\deps-vs2017-x86-installed\mpfr                                 ^
          -DMPFR_LIBRARY_DIR=%IFCOPENSHELL_ROOT%\deps-vs2017-x86-installed\mpfr                                 ^
          ..


License
-------

LGPL

Example
-------

in:

![](examples/rects.svg)

out:

![](examples/rects_output.svg)
