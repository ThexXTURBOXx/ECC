cd cocoalib
./configure --only-cocoalib
make lib doc -j4 LDFLAGS="-static"
make install -j4
cd ..
