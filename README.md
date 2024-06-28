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
using namespace O3DS;
<b>Modules:</b>
<ul>
  <li>#include "Open3DScan/general.hpp" - all moduls</li>
  <li>#include "Open3DScan/streams.hpp" - stream interface and std wrapper</li>
  <li>#include "Open3DScan/image.hpp" - image implement</li>
  <li>#include "Open3DScan/img_formats.hpp" - image formats interface and <b>BMP</b> implement</li>
  <li>#include "Open3DScan/linear_math.hpp" - matrix and vector</li>
  <li>#include "Open3DScan/numerics.hpp" - numerics methods</li>
  <li>#include "Open3DScan/formats_3D.hpp" - 3D-model format interface and <b>PLY</b> implement</li>
  <li>#include "Open3DScan/image_processing.hpp" - blur, canny, line detection, rectangle detection</li>
  <li>#include "Open3DScan/calculate_3D.hpp" - camera calibration and point cloud processing</li>
</ul>
To work with an image, use this class:
class img_t;<br>
To load and save images using the O3DS::img_t class not in the <b>BMP</b> format, you must independently implement algorithms for reading and writing from the O3DS::streams::istream_interface_t and O3DS::streams::ostream_interface_t stream, inheriting from the O3DS::img_formats::format_interface_t structure<br>
The same applies to 3D model formats. His format interface is O3DS::formats_3D::format_interface_t<br>
To calibrate cameras, create an instance of the O3DS::proc_3d::cameras_pair_t class. Its constructor accepts two images. They should show a hollow parallelepiped without rotation relative to both cameras. The geometric parameters of the parallelepiped must be known.<br>
The rest should be obvious, good luck!


