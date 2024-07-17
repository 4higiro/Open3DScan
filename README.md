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
// Read/write image interface for other image formats (has only BMP implimentation)
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
// Write point cloud interface for other model formats (has only PLY implimentation)
#include <Open3DScan/formats_3D.hpp>
```

The image processing module will work faster if the platform contains a parallel computing device, but there is also a simplified version for the CPU. To enable CPU calculation mode, switch the global flag to false

```cpp
// Disable parallel computing
O3DS::compute::gpu_mode = false;
```

To obtain a point cloud from a stereopair of images with line highlights, you need to:

```cpp
// Create image object
O3DS::img_t left(L"left camera image path"), right(L"right camera image path");
```

If you need to load pictures in any other format, then define the type of this format, inheriting from format interface. For example:
```cpp
// JPEG format type definition
class jpg_t : public O3DS::img_formats::format_interface_t;
```

Similarly, upload pictures of a hollow box, the far and near planes of which should be visible. Create a stereo pair object from these two images. This is necessary to calculate the focal length of the cameras used. It is worth determining the physical parameters of the cameras and box in advance:

```cpp
// Camera parametrs
struct camera_data_t final
{
	const img_t& source; // Image with box
	const proc_2d::rectangle_detector_t detector;
	math::vec2 pixel_size; // Physical pixel size
};
// Box parametrs
struct calibration_data_t final
{
	double width = 1.0;
	double height = 1.0;
	double depth = 1.0;
};

// Cameras pair initialisation
O3DS::proc_3d::cameras_pair_t(const camera_data_t& cam_a, const camera_data_t& cam_b, const calibration_data_t& data) noexcept;
```

Then upload pictures of highlights containing at least three square marks that the square detector can recognize. This is necessary to determine the rotations and relative positions of the cameras:

```cpp
// Incapsulation all camera calibration parametrs
struct cameras_proc_data_t final
{
	math::vec3 cama_euler, camb_euler;
	math::vec3 cama_position, camb_position;
	double cama_h = 0.0, camb_h = 0.0;
};

// Markers detection:
O3DS::proc_2d::rectangle_detector_t:get_line(const img_t& img, std::list<math::vec2>& points) const noexcept;

// Calculation parametrs in:
O3DS::proc_3d::cameras_pair_t::get_camera_rotation(proc_2d::rectangle_detector_t::rect_t source, math::vec3& target,
  camera_selection_t pick = camera_selection_t::A) const noexcept;

O3DS::proc_3d::camers_pair_t::get_cameras_position(std::vector<proc_2d::rectangle_detector_t::rect_t> rects_a,
  std::vector<proc_2d::rectangle_detector_t::rect_t> rects_b,
  math::vec3 cama_euler, math::vec3 camb_euler, math::vec3& cama_position, math::vec3& camb_position) noexcept;
```

Then create a handler object and build a point cloud. To do this, define the line detector object in advance:

```cpp
// Processor initialisation
O3DS::proc_3d::processor_t(const cameras_proc_data_t& cameras, const proc_2d::line_detector_t& detector,
  double dist_threshold = math::positive_inf) noexcept;

// Point clud processing
O3DS::proc_3d::processor_t::expand_point_cloud(const img_t& a, const img_t& b, std::vector<math::vec3> point_cloud) const noexcept;
```
