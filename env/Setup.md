# Setup

## On a fresh AWS Cloud9 Instance

1. Create a new AWS Cloud9 instance. Choose a configuration that gives you ~64GB of RAM (e.g. m4.4xlarge). Use a Amazon Linux 2 machine.

2. After starting up the instance, git clone this repository using HTTPS and authenticate.
On the command line, use:

`git clone https://github.com/Arcadia-Science/glial-origins.git`

3. Resize the local volume to enable more storage.
Set the number of GB of storage you want to have as an integer (e.g. 100GB below).
The resize script is from [this tutorial](https://docs.aws.amazon.com/cloud9/latest/user-guide/move-environment.html#move-environment-resize).
You probably want at least 100GB for a base analysis.
Up to 1000GB is not unreasonable if you're running a variety of analyses.

`bash env/resize.sh 100`

4. Download and install Miniconda.
On the command line, use:

`curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

`bash Miniconda3-latest-Linux-x86_64.sh`

> Accept all defaults.

5. Start Conda.
On the command line, use:

`source ~/.bashrc`

6. Set Conda base channels.
On the command line, use:

`conda config --add channels defaults`

`conda config --add channels bioconda`

`conda config --add channels conda-forge`


7. Install Mamba.
On the command line, use:

`conda install mamba`

Choose `y` when prompted.

8. Make sure you're in the `glial-origins/` directory.
For example, you might use:

`cd glial-origins/`

9. Create base environment.

`mamba env create -f env/glial_origins.yml`

## On an existing Cloud9 Instance

1. After starting up the instance, git clone this repository using HTTPS and authenticate.
On the command line, use:

`git clone https://github.com/Arcadia-Science/glial-origins.git`

2. If Conda or Mamba are not already installed, install them following steps 4-7 above.

3. Create the base environment using steps 8-9 above.
