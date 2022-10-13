cd cocoalib
./configure
make lib -j4 LDFLAGS="-static"
make install -j4
cd ..
