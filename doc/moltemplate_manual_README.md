Where can I find the moltemplate manual?
===========

### The PDF file for the moltemplate referance manual can be downloaded here:

**[https://www.moltemplate.org/doc/moltemplate_manual.pdf](https://www.moltemplate.org/doc/moltemplate_manual.pdf)**

Please let me know if the web site is down or the PDF file is out of date.


***Why?***  
I became concerned about the github download size.
The moltemplate manual is a large binary file which changes rapidly.
Every time a minor change was made to the manual,
the git repository would grow by about 700kB in size.



#### *Instructions for building the manual*

*Although there is no reason for most users to do this,
it is possible to construct the PDF file using the 
files in the *moltemplate_manual_src* subdirectory.
If you have *pdflatex* installed, you can use these commands to
create the "moltemplate_manual.pdf" file:*
```
cd moltemplate_manual_src/
pdflatex moltemplate_manual
pdflatex moltemplate_manual
bibtex moltemplate_manual
bibtex moltemplate_manual
pdflatex moltemplate_manual
pdflatex moltemplate_manual
mv moltemplate_manual.pdf ../
cd ../
```
