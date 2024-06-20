# this scripts builds the EAPM plugin into .hp format

# Enter the Plugin directory
cd EAPM

echo "Did you remove all deps from the deps folder EXCEPT bioprospecting? Y/n"

read isDeleted

if [ $isDeleted == "Y" -o $isDeleted == "y" ]
then

    # Read the OS_NAME (macOS or Linux)
    OS_NAME=$(uname)

    if [ $OS_NAME == "Darwin" ]
    then
        # Build the plugin for macOS
        echo "Building for macOS"
        OS_NAME="macOS"
    elif [ $OS_NAME == "Linux" ]
    then
        # Build the plugin for Linux
        echo "Building for Linux"
        OS_NAME="Linux"
    fi

    # Build the pages
    echo "Building Pages"
    npm run build:pages

    # Check that the EAPM/Pages folder exists
    if [ ! -d "Pages" ]
    then
        echo "Pages folder not found. Did the build fail?"
        exit 1
    fi

    # Remove previous build
    rm EAPM-$OS_NAME.hp

    # Zip the files
    zip -r EAPM-$OS_NAME.hp *

    echo "Plugin build finished"
else
    echo "Please remove them"
fi
