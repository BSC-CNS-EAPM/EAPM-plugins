# this scripts builds the EAPM plugin into .hp format

# Enter the Plugin directory
cd EAPM

# Zip the files
zip -r EAPM.hp *

# Move the zip file to the root directory
mv EAPM.hp ../