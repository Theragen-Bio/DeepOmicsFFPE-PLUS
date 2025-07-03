#!/bin/bash

miniconda_install_path="$HOME/Miniconda3-latest-Linux-x86_64.sh"

# 파일이 존재하지 않으면 wget 명령어 실행
if [ ! -f "$miniconda_install_path" ]; then
	echo "File does not exist. Downloading..."
	wget -O "$miniconda_install_path" "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
else
	echo "File already exists. Skipping download."
fi

if [ ! -d "$HOME/miniconda3_doffpe" ]; then
	echo "Miniconda not found. Installing..."
	/bin/bash "$miniconda_install_path" -p "$HOME/miniconda3_doffpe" -b
else
	echo "Miniconda path already exists. Skipping install."
fi

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('$HOME/miniconda3_doffpe/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
	eval "$__conda_setup"
else
	if [ -f "$HOME/miniconda3_doffpe/etc/profile.d/conda.sh" ]; then
		. "$HOME/miniconda3_doffpe/etc/profile.d/conda.sh"
	else
		export PATH="$HOME/miniconda3_doffpe/bin:$PATH"
	fi
fi
unset __conda_setup
# <<< conda initialize <<<

#conda info
if [ ! -d "$HOME/miniconda3_doffpe/envs/doffpe" ]; then
		echo "Creating conda environment doffpe..."
	conda create -y -n doffpe python=3.9.17
		conda activate doffpe
	conda install -y pysam==0.21.0 -c conda-forge -c bioconda
	pip install pandas==2.2.3
	pip install tqdm==4.66.1
	pip install requests pyYaml
	#pip install numpy==1.26.4
	#pip install einops==0.7.0
	#pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu102
	conda install -y -c bioconda htslib
	#pip install biopython==1.81
	conda install -y -c bioconda -c conda-forge bcftools=1.15.1
else
	echo "Conda environment doffpe already exists. Skipping creation."
fi

conda activate doffpe
conda info

PARENT_DIR="$(dirname "$(pwd)")"
chmod +x "$PARENT_DIR/doffpe"

# Check if the path is already in ~/.bashrc
if ! grep -q "$PARENT_DIR" ~/.bashrc; then
	echo "" >> ~/.bashrc
	echo "# Add DEEPOMICS FFPE PLUS tool path" >> ~/.bashrc
	echo "export PATH=\"\$PATH:$PARENT_DIR\"" >> ~/.bashrc
	echo "✅ 'doffpe' command has been registered. Please open a new terminal or run 'source ~/.bashrc' to apply the change."
else
	echo "ℹ️ The path is already registered in your .bashrc."
fi
