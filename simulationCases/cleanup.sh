#!/bin/bash

# cleanup.sh - Safely remove directories with confirmation
# WARNING: This script permanently deletes directories and all contents!
# Usage: ./cleanup.sh <directory_name>

# Check if argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No directory specified"
    echo "Usage: $0 <directory_name>"
    exit 1
fi

# Check if directory exists
if [ ! -d "$1" ]; then
    echo "Error: Directory '$1' does not exist"
    exit 1
fi

# Get absolute path for clarity
FULLPATH=$(cd "$1" 2>/dev/null && pwd)
if [ -z "$FULLPATH" ]; then
    echo "Error: Cannot access directory '$1'"
    exit 1
fi

# Show directory size and confirm deletion
echo "WARNING: This will permanently delete the following directory and ALL its contents:"
echo "  $FULLPATH"
echo ""
echo "Directory information:"
du -sh "$1" 2>/dev/null || echo "  Size: Unable to determine"
echo "  Number of files: $(find "$1" -type f 2>/dev/null | wc -l)"
echo ""
read -p "Are you sure you want to delete this directory? Type 'yes' to confirm: " CONFIRM

if [ "$CONFIRM" = "yes" ]; then
    echo "Deleting directory..."
    rm -rf "$1"
    if [ $? -eq 0 ]; then
        echo "Directory successfully deleted."
    else
        echo "Error: Failed to delete directory."
        exit 1
    fi
else
    echo "Deletion cancelled."
    exit 0
fi
