#!/bin/bash

# Run Doxygen to generate the documentation
doxygen Doxyfile

# Check if the symbolic link already exists and remove it if it does
if [ -L "index.html" ]; then
    rm index.html
fi

# Create the symbolic link to the index.html in the docs/html folder
ln -s docs/html/index.html index.html

echo "Documentation generated and symbolic link created in the root directory."

