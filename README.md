<h1>Open3DScan - is open source C++ framework (C++ 20 current)</h1>
<h2>Platforms: all</h2>
<h2>Libraries:</h2>
<ul>
  <li>OpenCL standart 1.1</li>
  <li>boost::compute 1.85.0</li>
  <li>OpenMP 2.0</li>
</ul>
<hr>
All tools are in the O3DS namespace<br>

```cpp
using namespace O3DS;
```

<h2>Moduls</h2>

```cpp
// All modules
#include <Open3DScan/general.hpp>
```

```cpp
// Base operations with image (scaling, rotation, crop and more...)
#include <Open3DScan/image.hpp>
```

```cpp
// Read/write image interface for other image formats (head only BMP implimentation)
#include <Open3DScan/image_formats.hpp>
```

```cpp
// Methods of image processing: convolution, blur, canny, line detection, marker detection
#include <Open3DScan/image_processing.hpp>
```

```cpp
// Math and numerics methods
#include <Open3DScan/numerics.hpp>
```

```cpp
// Matrix, vector and linear opeartion implementation
#include <Open3DScan/linear_math.hpp>
```

```cpp
// Common read/write interface + path processing
#include <Open3DScan/streams.hpp>
```

```cpp
// Camera calibration and point cloud building
#include <Open3DScan/calculate_3D.hpp>
```

```cpp
// Write point cloud interface for other model formats (head only PLY implimentation)
#include <Open3DScan/formats_3D.hpp>
```
