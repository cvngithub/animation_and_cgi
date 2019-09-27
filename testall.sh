files=(t1m1 t1m2 t1m3 t2m1 t2m2)
root=`pwd`

echo "STARTER_CODE TESTS"
for f in ${files[@]}; do
	cd $f
	rm -r test
	cp -r starter_code test
	cd test
	mkdir build
	cd build
	cmake -s -DCMAKE_CXX_FLAGS="" .. > /dev/null
	make --silent
	cd ..
	echo "$f STARTER CODE RESULTS"
	cp $root/run_tests.py .
	python3 run_tests.py $f
	echo "Press enter to continue..."
	read varname
	cd $root
	rm -rf $f/test
done

echo "SOLUTION TESTS"
for f in ${files[@]}; do
	cd $f
	rm -r test
	cp -r starter_code test
	cp solutions/* test/FOSSSim
	cd test
	mkdir build
	cd build
	cmake -s -DCMAKE_CXX_FLAGS="" .. > /dev/null
	make --silent
	cd ..
	cp $root/run_tests.py .
	echo "$f SOLUTION RESULTS"
	python3 run_tests.py $f
	echo "Press enter to continue..."
	read varname
	cd $root
	rm -rf $f/test
done
