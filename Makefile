
.PHONY: realigner
realigner:
	test -d venv3 || python3 -m venv venv3 --prompt "(realigner) "
	. ./venv3/bin/activate
	python3 setup.py build_ext --inplace

.PHONY: clean
clean:
	rm -r ./build
	rm -r ./__pycache__
	rm ./*.c
	rm ./*.so
	rm ./*.html

.PHONY: annotations
annotations:
	cython -a --3str aln.pyx
	cython -a --3str cig.pyx
