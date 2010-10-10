cd ..
./configure-build-install-default.sh
cd test-external-zip
rm out.* test1.exe dummy.txt.* test1.cxx.yaatk
g++ test1.cxx -I.. -L../_build-default/mdtk -lmdtk -lz -o test1.exe
./test1.exe
