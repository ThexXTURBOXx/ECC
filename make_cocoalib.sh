cd cocoalib
./configure --only-cocoalib
make lib -j4 LDFLAGS="-static"
make install -j4
cd ..
