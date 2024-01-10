# dune-pd-ana
A C++ analysis framework to estimate the performance of the DUNE Photon Detection System
Here you can find macro (file.cpp in folder Macro) useful to plot the
persistance of a set of waveform, compute the Signal-to-Noise ratio from
a calibration run, compute the FFTs of the electronic noise... In
Macro_description.txt you can find a brief overview of what the macros do.
"Header" contains the .hpp file where all the functions are delcared,
subdivided in categories. The flc.hpp file is the only one you should change to
adapt the codes to your dataset: it contains all the constant values needed for
your analysis.

To run a macro, you can use the command "root macro.cpp"


Federico Galizzi
