Generate Huge Fits Cube
=======================
This script generates an empty dummy 4 dimensional data cube in fits format.
After the initialisation this cube gets filled with fits image data. The
cube header gets updated from the first image in PATHLIST_STOKESI.
This script can be used to generate fits data cubes of sizes that exceeds the
machine's RAM (tested with 234 GB RAM and 335 GB cube data).

Please adjust the INPUT section in the script to your needs.

The data in directory `images` is test data and consists of Gaussian noise only.  

Developed at: IDIA (Institure for Data Intensive Astronomy), Cape Town, ZA  
Inspired by: https://github.com/idia-astro/image-generator  
Source: https://github.com/lh-astro/IDIA-generate-huge-fits-cube
E-Mail: astro[et]lennartheino[.]de
