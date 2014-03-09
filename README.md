Embedded Zerotree Wavelet Compression  
Tim Thirion  
Fall 2006

This is an implementation of the zerotree wavelet compression algorithm. It's
old and slow but the code is easy to follow. It should be useful for anyone
wanting to understand this algorithm.

See the data set chapelhill.txt for the correct file format.

To compile:
```
clang++ main.cpp -o ezw
```

To compress:
```
./ezw -c chapelhill.txt ch.ezw
```

To decompress:
```
./ezw -d ch.ezw reconstructed.txt
```
