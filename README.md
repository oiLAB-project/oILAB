# oILAB

Documentation
-------------
https://oilab-project.github.io/oILAB/


Requirements: Eigen 3.4 or later versions

To Build
--------

1) mkdir build
2) cd build
3) cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Release ..
4) make --build . --parallel

To Install
----------

After building, you can install oILAB system-wide or to a custom location:

### Custom installation location:
```bash
mkdir build
cd build
cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install ..
cmake --build . --parallel --target install
```

### Using oILAB in your CMake project:
After installation, you can use oILAB in your CMake project by adding:
```cmake
find_package(oILAB REQUIRED)
target_link_libraries(your_target PRIVATE oILAB::oILAB)
```

In your C++ code, include headers using:
```cpp
#include <oILAB/Lattices/LatticeModule.h>
#include <oILAB/Math/RationalMatrix.h>
```
