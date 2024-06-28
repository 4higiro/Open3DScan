#include "formats_3D.hpp"
#include "streams.hpp"
#include <string>

O3DS::formats_3d::format_interface_t::format_interface_t(std::vector<math::vec3>& cloud_to_save) noexcept
	: cloud(cloud_to_save) {}

O3DS::formats_3d::ply_t::ply_t(std::vector<math::vec3>& cloud_to_save) noexcept
 : O3DS::formats_3d::format_interface_t(cloud_to_save) {}

bool O3DS::formats_3d::ply_t::save(streams::ostream_interface_t& stream) const noexcept
{
	size_t vert_count = cloud.size();
	stream.write("ply\n", 0);
	stream.write("format ascii 1.0\n", 0);
	stream.write(("element vertex " + std::to_string(vert_count)).c_str(), 0);
	stream.write("\nproperty float x\n", 0);
	stream.write("property float y\n", 0);
	stream.write("property float z\n", 0);
	stream.write("element face 0\n", 0);
	stream.write("property list uchar int vertex_indices\n", 0);
	stream.write("end_header\n", 0);
	for (const auto& vert : cloud)
	{
		stream.write(std::to_string(vert["x"]).c_str(), 0);
		stream.write(" ", 0);
		stream.write(std::to_string(vert["y"]).c_str(), 0);
		stream.write(" ", 0);
		stream.write(std::to_string(vert["z"]).c_str(), 0);
		stream.write("\n", 0);
	}
	return true;
}
