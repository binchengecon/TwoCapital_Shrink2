# TwoCapital
Repo for R&amp;D model

## Prerequisite

Install PETSc, Python3.

## Intalling PETSc 

PETSc requires a Linux or Unix based OS. For Windows machine we therefore need to simulate such an OS. For MacOS machines unles they contain an M1 chip (M1 CPU) one can configure (install) the library directly following the instructions provided on the PETSc webpage. If the MacOS machine on which an installation is required contains a CPU in the M1 family please follow the instructions below. 

The basic idea is to create a virtual OS environment either by using the standard (and old) method of creating virtual machines or, alternatively, using docker containers. Virtual machines (VMs) tend to consume more resources compared to docker containers which are much more light-weight and, with the correct set up, are generally faster. Virtual machines are, however, overall more user friendly as the "correct set up" of docker containers requires some technical workarounds. 


### <a name="windows"></a>Installing on Windows

1. Download and install VirtualBox
2. Download an image file for the operating system to be simulated e.g. Ubuntu is a good choice and the latest version (20.04.03 as of Feb 8, 2022) is available [here](https://ubuntu.com/download/desktop).
3. Through virtual box create a virtual machine for a Linux OS. Make sure to leave enough space for the virtual machine i.e. when prompted assign 30GB rather than the default 8GB of disk space and assign at least 8GB (8192 MB) of RAM.
4. In terminal run: 

  `$ sudo apt-get install build-essential`

  this will install a `C++` compiler and other essentials

5. In terminal run: 

  `$ sudo apt-get install gfortran`

  this will install a `FORTRAN` compiler.

6. Install `Cython` - `Python` bindings for `C`. This can be done in 2 ways:
    - Installing the `Cython` [package](https://cython.readthedocs.io/en/stable/src/quickstart/install.html) directly with

    `$ pip install Cython` 
    - Installing `Anaconda` which contains `Cython` within the standard distribution. Installing `Anaconda` has advantages. It's a complete package management system which provide all necessary libraries and/or packages for scientific computing. `Anaconda` does not replace the default package manager `pip` but rather complements it and allows for environment management i.e. creating different workspaces containing different `Python` packages to enable easy testing, backwards compatibility checks, and resolution of package dependency issues. I user friendly installation tutorial on Ubuntu can be found [here](https://phoenixnap.com/kb/how-to-install-anaconda-ubuntu-18-04-or-20-04).

7. Follow the instructions on `Python` bindings [below](#python_bindings). 



### Installing on MacOS (M1 CPU)
The idea is to coerce docker containers to behave like VMs. This can be achieved through different VM management software such as Docker (Docker Destkop) or Vagrant:

For Docker:

1. Download and install Docker Desktop
2. Download an Ubuntu [image](https://hub.docker.com/_/ubuntu)
3. Create an Unbutu container with Docker
4. Follow steps 4 through 7 from the section ['Installing on Windows'](#windows).


For Vagrant:
1. Follow the instructions in this [tutorial](https://medium.com/nerd-for-tech/developing-on-apple-m1-silicon-with-virtual-environments-4f5f0765fd2f).


### <a name="python_bindings"></a> PETSc  with `Python` bindings

To run PETSc with `Python` we also need `Python` bindings provided by the `petsc4py` package. To install `petsc4py` and properly configure PETSc to work with `Python` follow these steps:

1. In the folder where PETSc download folder run the following configuration command: 

```
  $ ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack --with-shared-libraries`
  $ make all check
```

2. Install `petsc4py` following the official [documentation](https://github.com/erdc/petsc4py/blob/master/docs/source/install.rst).

The relevant part is: 

```
  $ export PETSC_DIR=/path/to/petsc
  $ export PETSC_ARCH=arch-linux2-c-opt
  $ pip install petsc4py
```

Make sure to change the environment variables `PETSC_DIR` and `PETSC_ARCH` to the ones obtained from the configuration command above (they are printed out on the terminal feed).

3. Configure again with: 

```
  $ ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack --with-shared-libraries --with-petsc4py
  $ make all check 
```








## Install the linear solver and other necessary packages


- `pip install -r requirements.txt`
- To use C implementation `pip install ./src/linearsystemcore`
- (Optional) installation of `SolveLinSys` package `pip install ./src/cppcore`

## Scripts and Model

A write-up is under the `./write-ups/write-up.pdf`

In py file, `linearsolver` stands for the linear solver to use: 
available solution:

- `pestsc` for PETSc + C implementation of coefficient matrix
- `petsc4py` for PETSc and numpy sparse matrix
- `eigen` for `SolveLinSys`



`./post-jump/post-jump-change.py` corresponds to section 1.2.1 in `write-up.pdf` with state variables log K, R, Y

`./tech4D/HJB-4d.py` corresponds to section 2 Pre jump HJB in `write-up.pdf`


