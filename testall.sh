regex=t[1-2]*
root=`pwd`

echo "STARTER_CODE TESTS"
for f in $regex; do
	cd $f
	cp -r starter_code test
	cp solutions/* test/FOSSSim
	cd test
	mkdir build
	cd build
	cmake .. && make
	cd ..
	python3 run_tests.py $f
	cd $root
	rm -rf $f/test
done

echo "SOLUTION TESTS"
for f in $regex; do
	cd $f
	cp -r starter_code test
	cd test
	mkdir build
	cd build
	cmake .. && make
	cd ..
	python3 run_tests.py $f
	cd $root
	rm -rf $f/test
done
