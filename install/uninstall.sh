#!/bin/bash

# === Step 1: Define installation-related paths ===
miniconda_install_path="$HOME/Miniconda3-latest-Linux-x86_64.sh"
dofp_miniconda_dir="$HOME/miniconda3_dofp"
dofp_cache_dir="$HOME/.dofp"

echo "Removing Miniconda installer and DOFP cache directory..."

# Remove Miniconda installer and local DOFP cache directory
rm -rf $miniconda_install_path $dofp_miniconda_dir $dofp_cache_dir

# === Step 2: Get parent directory path of current directory ===
PARENT_DIR="$(dirname "$(pwd)")"

# === Step 3: Backup current .bashrc ===
echo "Backing up .bashrc to prevent accidental loss..."
cp $HOME/.bashrc $HOME/.bashrc.backup.$(date +%Y%m%d_%H%M%S)

# === Step 4: Escape slashes for use in sed ===
ESCAPED_DIR=$(echo "$PARENT_DIR" | sed 's/\//\\\//g')

# === Step 5: Remove DOFP PATH entries from .bashrc ===
echo "Cleaning up .bashrc by removing DOFP PATH entries..."
sed -i '/# Add DOFP tool path/d' $HOME/.bashrc
sed -i "/export PATH=\\\"\$PATH:$ESCAPED_DIR\\\"/d" $HOME/.bashrc
echo "Uninstallation cleanup complete."
