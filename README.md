Moltemplate
===========

##  Description

Moltemplate is a cross-platform text-based molecule builder for LAMMPS. 

## Typical usage

    moltemplate.sh [-atomstyle style] [-pdb/-xyz coord_file] [-vmd] system.lt

## Web page

Documentation, examples, and supporting code can be downloaded at:

http://www.moltemplate.org

## Requirements

Moltemplate requires the Bourne-shell, and a recent version of python 
(2.7, 3.0 or higher), and can run on OS X, linux, or windows. (...if a 
suitable shell environment has been installed.  See below.)


## INSTALLATION INSTRUCTIONS

This directory should contain two folders:

    src/                          <-- contains all python and bash scripts
    examples/                     <-- simple examples built with moltemplate

Note: Popular force fields and molecule types are located here:

    src/moltemplate_force_fields/

The `moltemplate.sh` script and the python scripts that it invokes are 
located in the `src/` subdirectory.  You should update your PATH environment 
variable to include this directory.  

An alternate way to install moltemplate is to simply copy the entire contents
of moltemplate’s `src` subdirectory into a directory which is already in your
path, such as `/usr/local/bin`. (This will require admin privileges.)

### Installation example

Suppose the directory with this README.TXT file is located at ~/moltemplate.

If you use the bash shell, typically you would edit your 
`~/.profile`, `~/.bash_profile` or `~/.bashrc` files to contain the following line:

    export PATH="$PATH:$HOME/moltemplate/src"

If you use the tcsh shell, typically you would edit your 
`~/.login`, `~/.cshrc`, or `~/.tcshrc` files to contain the following lines:

    setenv PATH "$PATH:$HOME/moltemplate/src"

If you do not know what a `PATH` environment variable is and are curious, read:
    http://www.linfo.org/path_env_var.html
(I receive this question often.)

### WINDOWS installation suggestions

   You can install both moltemplate and LAMMPS in windows, but you will first need to install the BASH shell environment on your computer.  If you are using Windows 10 or later, try installing "Bash on Ubuntu on Windows":

https://msdn.microsoft.com/en-us/commandline/wsl/install_guide

   I have not tried this.  If that fails, try following the tutorial written by Yanqing Fu:

https://sourceforge.net/p/lammps/mailman/message/32599824/

   To use LAMMPS and moltemplate, you will also need to install (and learn how to use) a text editor.  (Word, Wordpad, and Notepad will probably not work.)  Popular text editors include: emacs, xemacs, vim, vi, jove, nano, gedit.  A longer list can be found here:
https://en.wikipedia.org/wiki/List_of_text_editors

## License

Moltemplate is available under the terms of the open-source 3-clause BSD 
license.  (See `LICENSE.txt`.)

