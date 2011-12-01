all: debug_static 

debug_shared: CONFIGURATION=debug_shared
debug_shared: libs
debug_static: CONFIGURATION=debug_static
debug_static: libs
release_shared: CONFIGURATION=release_shared
release_shared: libs
release_static: CONFIGURATION=release_static
release_static: libs

libs:
	cd online_structured_svm; make $(CONFIGURATION)
	cd blob; make $(CONFIGURATION)
	cd behavior; make $(CONFIGURATION)

clean:
	cd online_structured_svm; make clean
	cd blob; make clean
	cd behavior; make clean
	
