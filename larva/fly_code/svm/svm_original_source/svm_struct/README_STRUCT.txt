Readme for the SVM-light structure learning API
-----------------------------------------------
Thorsten Joachims, 03.07.2004

The API allows to implement different versions of the learning
algorithm for learning different kinds of structures. To adapt to a
new structure, one needs to modify the files

    svm_struct_api_types.h
    svm_struct_api.c

Both files already contain empty templates. The first file contains
the type definitions that need to be changed. PATTERN is the structure
for storing the x-part of an example (x,y), LABEL is the y-part. The
learned model will be stored in STRUCTMODEL. Finally,
STRUCT_LEARN_PARM can be used to store any parameters that you might
want to pass to the function.

The second file contains the function you need to implement. See the
documentation in the file for details.
