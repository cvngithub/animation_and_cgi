files=(t1m1 t1m2 t1m3 t2m1 t2m2)
root=`pwd`

echo "STARTER_CODE TESTS"
for f in ${files[@]}; do
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
