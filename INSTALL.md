# How to install FEniCS

## In the computer rooms in Jussieu:
To use git and fenics in the classroom:

At the command line open your `.bashrc` file in a text editor
```
subl $HOME/.bashrc
```
Then add these lines in the file  and save (modify with your student number and password)

```
git config --global  http.proxy http://studentnumber:password@134.157.103.2:3128
export https_proxy=http://studentnumber:password@134.157.103.2:3128
source /opt/fenics-2019/conda-fenics-2019.sh
export DIJITSO_CACHE_DIR=$HOME/.dijitsocache
mkdir -p $DIJITSO_CACHE_DIR 
```

## On your laptop

Instructions on how to install FEniCS are available https://fenicsproject.org/download/.
You also need a working `python3` installation with several additional packages. 

We suggest at first to try to use Anaconda, which as advanced package manager that can provide you without effort most of the software that you need. Install the Anaconda version based on *python3*.

### Anaconda install

This procedure will work on Mac and any Linux distribution.  Installation with Anaconda on windows should be tested. If you have Windows try it. If it does not work use a virtual machine.

First install Anaconda, following the instruction here  `https://docs.continuum.io/anaconda/install/`. 

Then run following commands in your terminal:

```
conda config --add channels conda-forge 
conda config --set channel_priority strict
conda create -n fenicsproject 
source activate fenicsproject
conda install fenics mshr
```

We suggest also to install other useful packages by executing the following. We use also `pip` which is the default python package manager:

```
conda install gmsh pip
pip3 install vtkplotter meshio jupyter spyder matplotlib scipy
```

If you need to install further package, just type `pip3 install packagename` or `conda install packagename`

**Important:** You need to have the conda `fenicsproject` environment activated when you install the package or use python. If the environment is activate you terminal prompt will include `(fenicsproject)`, like the following:
```
(fenicsproject) maurini@maurinis-MacBook-Pro:~$ 
```

## Cloud-based installations
You can run python programs jupyter notebooks and FEniCS on online servers. The basic service is free and can be a solution if all other installation systems fail.

* *advantages:* You do not need to install anything on your machine and you do note use the ressources of your machine

* *disadvantages:* It can be slow (especially Azure, colab seems better). You can use only jupyter notebooks, you share your date with google or microsoft, you do not a full control of the system, you need to be online with a good network connection.

### Google colab
* You need a google account
* You can use this example to start https://colab.research.google.com/drive/1yhw-E21eEa4-ej27abiJiHj064iRaetS#scrollTo=iF68ZCvDTvo5
* You can save the notebooks and your working environment on your google drive

### Microsoft Azure
Another possibility is using azure notebbok (you need a Microsoft account)
https://fenics-cmaurini.notebooks.azure.com/nb/notebooks/FEniCS%20demo.ipynb

## Other possibilities

### On ubuntu (preferred installation):
````
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install --no-install-recommends fenics
````

### Using docker virtual machines (suggested method if other fail).:
* See https://fenicsproject.org/download/ on the docker section and https://fenics-containers.readthedocs.io/en/latest/index.html
