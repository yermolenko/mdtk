cd ..
./configure-build-install-default.sh
cd test-in-memory-gzip
rm out.* test1.exe
g++ test1.cxx -I.. -L../_build-default/mdtk -lmdtk -lz -o test1.exe
./test1.exe
