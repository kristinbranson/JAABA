Code for Latent Structural SVM (June 4 2009)
--------------------------------------------
(implemented by Chun-Nam Yu based on the SVM-light code by Thorsten Joachims)


INSTALLATION
------------
1. Type 'make' in the directory latentssvm and make sure that the compilation is successful.


USING THE LATENT SVM^STRUCT API
-------------------------------
1. This implemention of latent structural SVM follows the modular design of SVM^struct [2]. To implement your own latent structural SVM application, you need to modify the two files svm_struct_latent_api_types.h and svm_struct_latent_api.c. 
2. An example implementation on the application of motif finding in [1] is available for download at http://www.cs.cornell.edu/~cnyu/latentssvm. 
3. You may also find the set of example applications on the SVM^struct website instructive when implementing your own latent structural SVM applications. The structural SVM application examples are available at http://svmlight.joachims.org/svm_struct.html.


CONTACT
-------
If you had any suggestions to the program or have bugs to report, you can email Chun-Nam Yu at cnyu@cs.cornell.edu.  


REFERENCES
----------
[1] C.-N. Yu and T. Joachims: Learning Structural SVMs with Latent Variables, ICML 2009 
[2] I. Tsochantaridis, T. Hofmann, T. Joachims, and Y. Altun: Support Vector Learning for Interdependent and Structured Output Spaces, ICML 2004
