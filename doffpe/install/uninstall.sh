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
#ESCAPED_DIR=$(echo "$PARENT_DIR" | sed 's/\//\\\//g')

# === Step 5: Remove DEEPOMICS FFPE PLUS PATH entries from .bashrc ===
#echo "Cleaning up .bashrc by removing doffpe PATH entries..."
#sed -i '/# Add DEEPOMICS FFPE PLUS tool path/d' $HOME/.bashrc
#sed -i "/export PATH=\\\"$ESCAPED_DIR:\$PATH\\\"/d" $HOME/.bashrc
#echo "Uninstallation cleanup complete."

# === Step 4: Get parent directory ===
PARENT_DIR="$(dirname "$(pwd)")"

# === Step 5: Escape slashes for use in sed ===
ESCAPED_DIR=$(echo "$PARENT_DIR" | sed 's/\//\\\//g')

# === Step 6: Remove DEEPOMICS FFPE PLUS PATH entries from .bashrc ===
echo "🧹 Cleaning up .bashrc by removing doffpe PATH entries..."
sed -i '/# Add DEEPOMICS FFPE PLUS tool path/d' "$HOME/.bashrc"
sed -i "/export PATH=\\\"$ESCAPED_DIR:\\\$PATH\\\"/d" "$HOME/.bashrc"

# === Step 7: Also remove from current shell PATH if present ===
if [[ ":$PATH:" == *":$PARENT_DIR:"* ]]; then
    export PATH=$(echo "$PATH" | sed -e "s|$PARENT_DIR:||" -e "s|:$PARENT_DIR||" -e "s|$PARENT_DIR||")
    echo "✅ Removed $PARENT_DIR from current shell PATH."
else
    echo "ℹ️ $PARENT_DIR not found in current shell PATH."
fi

echo "✅ Uninstallation cleanup complete."
