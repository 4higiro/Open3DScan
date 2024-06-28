#ifndef _OPEN_3D_SCAN_3D_FORMATS
#define _OPEN_3D_SCAN_3D_FOTMATS

#include "linear_math.hpp"
#include <vector>

namespace O3DS
{
	namespace streams { struct istream_interface_t; struct ostream_interface_t; }

	namespace formats_3d
	{
		struct format_interface_t
		{
			std::vector<math::vec3>& cloud;
			format_interface_t(std::vector<math::vec3>& cloud_to_save) noexcept;
			virtual bool save(streams::ostream_interface_t& stream) const noexcept = 0;
		};

		struct ply_t : format_interface_t
		{
			ply_t(std::vector<math::vec3>& cloud_to_save) noexcept;
			bool save(streams::ostream_interface_t& stream) const noexcept override;
		};
	}
}

#endif
