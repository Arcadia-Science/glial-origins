# Setup

## On a fresh AWS Cloud9 Instance

1. Create a new AWS Cloud9 instance.  
- Give it a useful name (e.g. `glial-origins-analysis`)
- Choose a configuration that gives you ~64GB of RAM (e.g. m4.4xlarge).  
- Use a Amazon Linux 2 machine.

2. After starting up the instance, git clone this repository using HTTPS and authenticate.  
On the command line, use:  
`git clone https://github.com/Arcadia-Science/glial-origins.git`. 
> Temporarily, you'll want to use the latest branch.  
> `git checkout das/orthofinder-dev`  
> In the next PR, you can do this install from `main`.

3. Move into the GitHub repo.  
`cd glial-origins/`

4. Resize the local volume to enable more storage.  
Set the number of GB of storage you want to have as an integer (e.g. 100GB below).  
The resize script is from [this tutorial](https://docs.aws.amazon.com/cloud9/latest/user-guide/move-environment.html#move-environment-resize).   
`bash env/resize.sh 100`
> You probably want at least 100GB for a base analysis.  
> Up to 1000GB is not unreasonable if you're running a variety of analyses. 

5. Download and install Miniconda.  
On the command line, use:  
`curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`  
`bash Miniconda3-latest-Linux-x86_64.sh`  
> Accept all defaults.

6. Start Conda.  
On the command line, use:  
`source ~/.bashrc`

7. Set Conda base channels.  
On the command line, use:  
`conda config --add channels defaults`  
`conda config --add channels bioconda`  
`conda config --add channels conda-forge`  

8. Install Mamba.  
On the command line, use:  
`conda install mamba`
> Choose `y` when prompted.

9. Create base environment.  
`mamba env create -f env/glial_origins.yml`  
> This will take some time, perhaps ~5-10min.

10. Activate the new conda environment.  
`conda activate glial_origins`

11. Run Jupyter.  
`jupyter lab --ip 0.0.0.0 --port 8888`
> This should output a bunch of lines of text.  
> You'll want to pay attention to the lines that say something like:
` To access the server, open this file in a browser:
        file:///home/ec2-user/.local/share/jupyter/runtime/jpserver-30829-open.html
    Or copy and paste one of these URLs:
        http://ip-172-31-22-55.us-west-1.compute.internal:8888/lab?token=dbaec68e69abc5a9e6ccc5bf413eb2a170c938c37ddb55c1
     or http://127.0.0.1:8888/lab?token=dbaec68e69abc5a9e6ccc5bf413eb2a170c938c37ddb55c1`
> You'll use the long string after `token=` to set the passphrase for the Jupyter lab server in step 13.

12. Go to the Cloud9 console and go to the EC2 instance.  
- From the Cloud9 "Your Environments" section, click on the name of your new environment (e.g. `glial-origins-analysis`).
- Under "EC2 Instance", click on "Go To Instance."
- Select the instance in the checkbox menu.
- On the bottom half of the page, click the "Security" tab.
- Click on the hyperlink under "Security Groups."
- Click the "Edit Inbound Rules" button at the top-right of the bottom table.
- Click the "Add Rule" button.
- Under "Port Range" type `8888`.
- Under "Source" choose `Anywhere - IPv4`.
- Click the orange "Save rules" button.

13. Login to the Jupyter server.
- From the EC2 Instance checkbox menu, select the instance.
- Copy the "Public IPv4 DNS" value.
- Paste into a new web browser tab and add `:8888` to the end of the address.
- Copy and paste the string after `token=` in step 11 into the `Token` box.
- Type a secure password into the `New Password` box.
- Click the `Log in and set new password` button.

14. Start working on analyses!

## On an existing Cloud9 Instance

1. After starting up the instance, git clone this repository using HTTPS and authenticate.  
On the command line, use:  
`git clone https://github.com/Arcadia-Science/glial-origins.git`  

2. If Conda or Mamba are not already installed, install them following steps 4-7 above.

3. Create and start the base environment using steps 8-10 above.

4. Run Jupyter lab and set up the server using steps 11-13 above. If you have already set security group rules, you may not need to do this again.

5. Start working on analyses!
