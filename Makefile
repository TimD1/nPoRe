.PHONY: npore
npore:
	python3 setup.py build_ext --inplace
	mv *.so src/

.PHONY: clean
clean:
	rm -rf ./build
	rm -rf src/__pycache__
	rm -f src/*.c
	rm -f src/*.so
	rm -f src/*.html

.PHONY: annotations
annotations:
	cython -a --3str src/aln.pyx
	cython -a --3str src/cig.pyx
	cython -a --3str src/bam.pyx
