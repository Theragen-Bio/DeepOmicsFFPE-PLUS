#!/bin/bash

# === Step 1: Define installation-related paths ===
miniconda_install_path="$HOME/Miniconda3-latest-Linux-x86_64.sh"
doffpe_miniconda_dir="$HOME/miniconda3_doffpe"
doffpe_cache_dir="$HOME/.doffpe"

echo "Removing Miniconda installer and doffpe cache directory..."

# Remove Miniconda installer and local doffpe cache directory
rm -rf $miniconda_install_path $doffpe_miniconda_dir $doffpe_cache_dir

# === Step 2: Get parent directory path of current directory ===
PARENT_DIR="$(dirname "$(pwd)")"

# === Step 3: Backup current .bashrc ===
echo "Backing up .bashrc to prevent accidental loss..."
cp $HOME/.bashrc $HOME/.bashrc.backup.$(date +%Y%m%d_%H%M%S)

# === Step 4: Escape slashes for use in sed ===
ESCAPED_DIR=$(echo "$PARENT_DIR" | sed 's/\//\\\//g')

# === Step 5: Remove DEEPOMICS FFPE PLUS PATH entries from .bashrc ===
echo "Cleaning up .bashrc by removing doffpe PATH entries..."
sed -i '/# Add DEEPOMICS FFPE PLUS tool path/d' $HOME/.bashrc
sed -i "/export PATH=\\\"\$PATH:$ESCAPED_DIR\\\"/d" $HOME/.bashrc
echo "Uninstallation cleanup complete."
