all: extensions

extensions:
	python setup.py build_ext --inplace
	rm -rf build/

clean:
	rm -f *.so
	rm -f *.pyc
	rm -f *.c

