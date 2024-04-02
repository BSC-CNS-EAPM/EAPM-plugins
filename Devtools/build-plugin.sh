# this scripts builds the EAPM plugin into .hp format

# Enter the Plugin directory
cd EAPM

echo "Did you remove all deps from the deps folder EXCEPT bioprospecting? Y/n"

read isDeleted

if [ $isDeleted == "Y" -o $isDeleted == "y" ]
then
    # Remove previous build
    rm EAPM.hp

    # Zip the files
    zip -r EAPM.hp *

    # Move the zip file to the root directory
    mv EAPM.hp ../

    echo "Plugin build finished"
else
    echo "Please remove them"
fi
